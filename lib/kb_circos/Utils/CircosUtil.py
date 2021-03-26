import errno
import json
import os
import subprocess
import sys
import time
import uuid
import zipfile
import copy
import shutil
import math
import pandas as pd

from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.MetagenomeUtilsClient import MetagenomeUtils
from installed_clients.ReadsUtilsClient import ReadsUtils

from random import seed
from random import randint
# seed random number generator
seed(1)

from shutil import copyfile


def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))


class CircosUtil:
    CIRCOS_BASE_PATH = '/Circos'
    CIRCOS_RESULT_DIRECTORY = 'circos_output_dir'
    MAPPING_THREADS = 16
    BBMAP_MEM = '30g'

    def __init__(self, config):
        self.callback_url = config['SDK_CALLBACK_URL']
        self.scratch = config['scratch']
        self.shock_url = config['shock-url']
        self.ws_url = config['workspace-url']
        self.dfu = DataFileUtil(self.callback_url)
        self.ru = ReadsUtils(self.callback_url)
        self.au = AssemblyUtil(self.callback_url)
        self.mgu = MetagenomeUtils(self.callback_url)

    def _validate_run_circos_params(self, task_params):
        """
        _validate_run_circos_params:
                validates params passed to run_circos method
        """
        log('Start validating run_circos params')

        # check for required parameters
        for p in ['assembly_ref', 'reads_file', 'workspace_name']:
            if p not in task_params:
                raise ValueError('"{}" parameter is required, but missing'.format(p))

    def _mkdir_p(self, path):
        """
        _mkdir_p: make directory for given path
        """
        if not path:
            return
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise

    def _run_command(self, command):
        """
        _run_command: run command and print result
        """
        os.chdir(self.scratch)
        log('Start executing command:\n{}'.format(command))
        log('Command is running from:\n{}'.format(self.scratch))
        pipe = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output, stderr = pipe.communicate()
        exitCode = pipe.returncode

        if (exitCode == 0):
            log('Executed command:\n{}\n'.format(command) +
                'Exit Code: {}\n'.format(exitCode))
        else:
            error_msg = 'Error running command:\n{}\n'.format(command)
            error_msg += 'Exit Code: {}\nOutput:\n{}\nStderr:\n{}'.format(exitCode, output, stderr)
            raise ValueError(error_msg)
            sys.exit(1)
        return (output, stderr)

    # this function has been customized to return read_type variable (interleaved vs single-end library)
    def stage_reads_file(self, reads_file):
        """
        stage_reads_file: download fastq file associated to reads to scratch area
                          and return result_file_path
        """

        log('Processing reads object list: {}'.format(reads_file))

        result_file_path = []
        read_type = []

        reads_file_check = isinstance(reads_file, list)
        if reads_file_check :
            log("Input reads_file detected as list. Great.")
        else:
            log("Input reads_file not a list, converting.")
            reads_file = [reads_file]


        # getting from workspace and writing to scratch. The 'reads' dictionary now has file paths to scratch.
        reads = self.ru.download_reads({'read_libraries': reads_file, 'interleaved': None})['files']

        # reads_file is the list of file paths on workspace? (i.e. 12804/1/1).
        # "reads" is the hash of hashes where key is "12804/1/1" or in this case, read_obj and
        # "files" is the secondary key. The tertiary keys are "fwd" and "rev", as well as others.
        for read_obj in reads_file:
            files = reads[read_obj]['files']    # 'files' is dictionary where 'fwd' is key of file path on scratch.
            result_file_path.append(files['fwd'])
            read_type.append(files['type'])
            if 'rev' in files and files['rev'] is not None:
                result_file_path.append(files['rev'])

        return result_file_path, read_type

    def _get_contig_file(self, assembly_ref):
        """
        _get_contig_file: get contig file from GenomeAssembly object
        """
        contig_file = self.au.get_assembly_as_fasta({'ref': assembly_ref}).get('path')

        sys.stdout.flush()
        contig_file = self.dfu.unpack_file({'file_path': contig_file})['file_path']

        return contig_file
    #
    def retrieve_assembly(self, task_params):
        if os.path.exists(task_params['contig_file_path']):
            assembly = task_params['contig_file_path']
            print("FOUND ASSEMBLY ON LOCAL SCRATCH")
        else:
            # we are on njsw so lets copy it over to scratch
            assembly = self._get_contig_file(task_params['assembly_ref'])
        return assembly

    def deinterlace_raw_reads(self, fastq):
        fastq_forward = fastq.split('.fastq')[0] + "_forward.fastq"
        fastq_reverse = fastq.split('.fastq')[0] + "_reverse.fastq"
        command = 'deinterleave_fastq.sh < {} {} {}'.format(fastq, fastq_forward, fastq_reverse)
        try:
            self._run_command(command)
        except:
            raise Exception("Cannot deinterlace fastq file!")
        return (fastq_forward, fastq_reverse)

    def run_read_mapping_interleaved_pairs_mode(self, task_params, assembly, fastq, sam):
        read_mapping_tool = task_params['read_mapping_tool']
        log("running {} mapping in interleaved mode.".format(read_mapping_tool))
        random_seed_int = randint(0, 999999999)
        log("randomly selected seed (integer) used for read mapping is: {}".format(random_seed_int))
        if task_params['read_mapping_tool'] == 'bbmap':
            log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
            command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
            command += 'threads={} '.format(self.MAPPING_THREADS)
            command += 'ref={} '.format(assembly)
            command += 'in={} '.format(fastq)
            command += 'out={} '.format(sam)
            command += 'fast interleaved=true mappedonly nodisk overwrite'
        elif task_params['read_mapping_tool'] == 'bowtie2_default':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            bt2index = os.path.basename(assembly) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 -x {} '.format(bt2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'bowtie2_very_sensitive':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            bt2index = os.path.basename(assembly) + '.bt2'
            command = 'bowtie2-build -f {} '.format(assembly)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} && '.format(bt2index)
            command += 'bowtie2 --very-sensitive -x {} '.format(bt2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '--threads {} '.format(self.MAPPING_THREADS)
            command += '-S {}'.format(sam)
        elif task_params['read_mapping_tool'] == 'minimap2':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            command = 'minimap2 -ax sr -t {} '.format(self.MAPPING_THREADS)
            command += '--seed {} '.format(random_seed_int)
            command += '{} '.format(assembly)
            command += '{} '.format(fastq_forward)
            command += '{} > '.format(fastq_reverse)
            command += '{}'.format(sam)
        elif task_params['read_mapping_tool'] == 'hisat2':
            (fastq_forward, fastq_reverse) = self.deinterlace_raw_reads(fastq)
            ht2index = os.path.basename(assembly) + '.ht2'
            command = 'hisat2-build {} '.format(assembly)
            command += '{} && '.format(ht2index)
            command += 'hisat2 -x {} '.format(ht2index)
            command += '-1 {} '.format(fastq_forward)
            command += '-2 {} '.format(fastq_reverse)
            command += '-S {} '.format(sam)
            command += '--seed {} '.format(random_seed_int)
            command += '--threads {}'.format(self.MAPPING_THREADS)
        log('running alignment command: {}'.format(command))
        out, err = self._run_command(command)
    #
    # def run_read_mapping_unpaired_mode(self, task_params, assembly, fastq, sam):
    #     read_mapping_tool = task_params['read_mapping_tool']
    #     log("running {} mapping in single-end (unpaired) mode.".format(read_mapping_tool))
    #     random_seed_int = randint(0, 999999999)
    #     log("randomly selected seed (integer) used for read mapping is: {}".format(random_seed_int))
    #     if task_params['read_mapping_tool'] == 'bbmap':
    #         log("Warning: bbmap does not support setting random seeds, so results are not reproducible.")
    #         command = 'bbmap.sh -Xmx{} '.format(self.BBMAP_MEM)
    #         command += 'threads={} '.format(self.MAPPING_THREADS)
    #         command += 'ref={} '.format(assembly)
    #         command += 'in={} '.format(fastq)
    #         command += 'out={} '.format(sam)
    #         command += 'fast interleaved=false mappedonly nodisk overwrite'
    #         # BBMap is deterministic without the deterministic flag if using single-ended reads
    #     elif task_params['read_mapping_tool'] == 'bowtie2_default':
    #         bt2index = os.path.basename(assembly) + '.bt2'
    #         command = 'bowtie2-build -f {} '.format(assembly)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '{} && '.format(bt2index)
    #         command += 'bowtie2 -x {} '.format(bt2index)
    #         command += '-U {} '.format(fastq)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '-S {}'.format(sam)
    #     elif task_params['read_mapping_tool'] == 'bowtie2_very_sensitive':
    #         bt2index = os.path.basename(assembly) + '.bt2'
    #         command = 'bowtie2-build -f {} '.format(assembly)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '{} && '.format(bt2index)
    #         command += 'bowtie2 --very-sensitive -x {} '.format(bt2index)
    #         command += '-U {} '.format(fastq)
    #         command += '--threads {} '.format(self.MAPPING_THREADS)
    #         command += '-S {}'.format(sam)
    #     elif task_params['read_mapping_tool'] == 'minimap2':
    #         command = 'minimap2 -ax sr -t {} '.format(self.MAPPING_THREADS)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '{} '.format(assembly)
    #         command += '{} > '.format(fastq)
    #         command += '{}'.format(sam)
    #     elif task_params['read_mapping_tool'] == 'hisat2':
    #         ht2index = os.path.basename(assembly) + '.ht2'
    #         command = 'hisat2-build {} '.format(assembly)
    #         command += '{} && '.format(ht2index)
    #         command += 'hisat2 -x {} '.format(ht2index)
    #         command += '-U {} '.format(fastq)
    #         command += '-S {} '.format(sam)
    #         command += '--seed {} '.format(random_seed_int)
    #         command += '--threads {}'.format(self.MAPPING_THREADS)
    #     log('running alignment command: {}'.format(command))
    #     out, err = self._run_command(command)

    def convert_sam_to_sorted_and_indexed_bam(self, sam):
        # create bam files from sam files
        sorted_bam = os.path.abspath(sam).split('.sam')[0] + "_sorted.bam"

        command = 'samtools view -F 0x04 -uS {} | '.format(sam)
        command += 'samtools sort - -o {}'.format(sorted_bam)

        log('running samtools command to generate sorted bam: {}'.format(command))
        self._run_command(command)

        # verify we got bams
        if not os.path.exists(sorted_bam):
            log('Failed to find bam file\n{}'.format(sorted_bam))
            sys.exit(1)
        elif(os.stat(sorted_bam).st_size == 0):
            log('Bam file is empty\n{}'.format(sorted_bam))
            sys.exit(1)

        # index the bam file
        command = 'samtools index {}'.format(sorted_bam)

        log('running samtools command to index sorted bam: {}'.format(command))
        self._run_command(command)

        return sorted_bam

    def generate_alignment_bams(self, task_params, assembly):
        """
            This function runs the selected read mapper and creates the
            sorted and indexed bam files from sam files using samtools.
        """

        reads_file = task_params['reads_file']

        (read_scratch_path, read_type) = self.stage_reads_file(reads_file)

        # sorted_bam_file_list = []

        # list of reads files, can be 1 or more. assuming reads are either type unpaired or interleaved
        # will not handle unpaired forward and reverse reads input as seperate (non-interleaved) files

        for i in range(len(read_scratch_path)):
            fastq = read_scratch_path[i]
            fastq_type = read_type[i]

            sam = os.path.basename(fastq).split('.fastq')[0] + ".sam"
            #sam = os.path.join(self.CIRCOS_RESULT_DIRECTORY, sam)

            if fastq_type == 'interleaved':  # make sure working - needs tests
                log("Running interleaved read mapping mode")
                self.run_read_mapping_interleaved_pairs_mode(task_params, assembly, fastq, sam)
            else:  # running read mapping in single-end mode
                log("Running unpaired read mapping mode")
                self.run_read_mapping_unpaired_mode(task_params, assembly, fastq, sam)

            sorted_bam = self.convert_sam_to_sorted_and_indexed_bam(sam)

        return sorted_bam

    def uppercase_fastq_file(self, reads_file):
        output_fastq = reads_file.rsplit('.', 1)[0] + "_uppercase.fastq"
        command = 'seqkit seq -u '
        command += '{} > '.format(reads_file)
        command += '{}'.format(output_fastq)
        log('uppercase_fastq_file: {}'.format(command))
        self._run_command(command)
        return output_fastq

    def clean_input_fasta(self, assembly):
        assembly_clean = assembly.rsplit('.', 1)[0] + "_clean.fasta"
        command = 'cut -d\' \' -f1 {} > {}'.format(assembly, assembly_clean)
        log('clean_input_fasta: {}'.format(command))
        self._run_command(command)
        return assembly_clean

    def sort_fasta_by_length(self, input_fasta):
        output_fasta_sorted = input_fasta.rsplit('.', 1)[0] + "_sorted.fasta"
        command = 'seqkit sort '
        command += '-l {} '.format(input_fasta)
        command += '-r '.format(input_fasta)
        command += '> {} '.format(output_fasta_sorted)
        log('sort_fasta_by_length: {}'.format(command))
        self._run_command(command)
        return output_fasta_sorted

    def extract_mapping_tracks_from_bam(self, sorted_bam):
        output_sam = 'final.mapped.sam'
        command = 'bedtools genomecov -ibam '
        command += '{} '.format(sorted_bam)
        command += '-bg > circos_mapping_tracks.txt'
        log('extract_mapping_tracks_from_bam: {}'.format(command))
        self._run_command(command)
        file1 = open(os.path.abspath("circos_mapping_tracks.txt"), 'r')
        pandas_df = pd.read_table(file1, header=None)
        max_cov = round(pandas_df[3].max(),1)
        min_cov = round(pandas_df[3].min(),1)
        std_cov = round(pandas_df[3].std(),1)
        mean_cov = round(pandas_df[3].mean(),1)
        return max_cov, min_cov, std_cov, mean_cov

    def make_circos_karyotype_file(self, assembly_clean_sorted):
        from Bio import SeqIO
        count = 1
        path_to_circos_karyotype_file = "circos_karyotype.txt"
        with open(path_to_circos_karyotype_file, 'a') as f:
            for record in SeqIO.parse(assembly_clean_sorted, "fasta"):
                f.write("chr - {} {} 0 {} {}\n".format(record.id, count, len(record), record.id))
                count += 1
        f.close()

    def count_num_contigs(self, assembly_clean_sorted):
        num_contigs = 0
        file1 = open(assembly_clean_sorted, 'r')
        lines = file1.readlines()
        for line in lines:
            if line.startswith('>'):
                num_contigs = num_contigs + 1
        return num_contigs


    def prep_circos_axis(self, max_cov):
        if max_cov < 30:
            max_cov = 30
        command = 'sed -i "s/^max.*/max   = {}/" /kb/module/lib/kb_circos/circos/circos.conf '.format(max_cov)
        log('prep_circos_axis: {}'.format(command))
        self._run_command(command)

    def draw_circos_plot(self):
        command = 'circos -conf '
        command += '/kb/module/lib/kb_circos/circos/circos.conf'
        log('draw_circos_plot: {}'.format(command))
        self._run_command(command)

    def make_circos_plot(self, task_params, reads_file, output_circos_assembly):
        output_fastq = self.uppercase_fastq_file(reads_file)
        output_circos_assembly_clean = self.clean_input_fasta(output_circos_assembly)
        output_circos_assembly_clean_sorted = self.sort_fasta_by_length(output_circos_assembly_clean)
        sam = os.path.basename(output_fastq).split('.fastq')[0] + ".sam"
        self.run_read_mapping_interleaved_pairs_mode(task_params, output_circos_assembly_clean_sorted, output_fastq, sam)
        sorted_bam = self.convert_sam_to_sorted_and_indexed_bam(sam)
        max_cov, min_cov, std_cov, mean_cov = self.extract_mapping_tracks_from_bam(sorted_bam)
        self.make_circos_karyotype_file(output_circos_assembly_clean_sorted)
        num_contigs = self.count_num_contigs(output_circos_assembly_clean_sorted)
        self.prep_circos_axis(max_cov)
        self.draw_circos_plot()
        return output_circos_assembly_clean_sorted, max_cov, min_cov, std_cov, mean_cov, num_contigs

    def move_circos_output_files_to_output_dir(self):
        dest = os.path.abspath(self.CIRCOS_RESULT_DIRECTORY)
        files = os.listdir(os.path.abspath(self.scratch))
        for f in files:
            if (f.startswith("circos") or f.endswith("sorted.fasta") or f.endswith("clean.fasta")):
                shutil.move(f, dest)

    def generate_circos_command(self, task_params):
        """
        generate_command: circos
        """

        assembly_ref = task_params['contig_file_path']
        reads_file = task_params['reads_list_file'][0]

        log("\n\nRunning generate_circos_command")

        assembly_clean_sorted, max_cov, min_cov, std_cov, mean_cov, num_contigs = self.make_circos_plot(task_params, reads_file, assembly_ref)

        self.move_circos_output_files_to_output_dir()

        return assembly_clean_sorted, max_cov, min_cov, std_cov, mean_cov, num_contigs



    def generate_output_file_list(self, result_directory):
        """
        generate_output_file_list: zip result files and generate file_links for report
        """
        log('Start packing result files')
        output_files = list()

        output_directory = os.path.join(self.scratch, str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file = os.path.join(output_directory, 'circos_result.zip')

        with zipfile.ZipFile(result_file, 'w',
                             zipfile.ZIP_DEFLATED,
                             allowZip64=True) as zip_file:

            for dirname, subdirs, files in os.walk(result_directory):
#                for file in files:
#                    zip_file.write(os.path.join(dirname, file), file)
                if (dirname.endswith(self.CIRCOS_RESULT_DIRECTORY)):
                    baseDir = os.path.basename(dirname)
                    for file in files:
                        full = os.path.join(dirname, file)
                        zip_file.write(full, os.path.join(baseDir, file))

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'Files generated by Circos App'})

        return output_files



    def generate_html_report(self, assembly_ref, assembly_stats):
        """
        generate_html_report: generate html summary report
        """

        log('Start generating html report')
        #html_report = list()

        output_directory = os.path.join(self.scratch, 'html_dir_' + str(uuid.uuid4()))
        self._mkdir_p(output_directory)
        result_file_path = os.path.join(output_directory, 'report.html')

        # get summary data from existing assembly object and bins_objects
        Summary_Table_Content = ''
        Overview_Content = ''

        # generate overview content

        # get png
        png_filename_l = [f for f in os.listdir(self.CIRCOS_RESULT_DIRECTORY) if f.endswith('.png')]

        # Example
        # Overview_Content += '<p>Input contigs: {}</p>'.format(input_contig_count)

        Overview_Content += '<p>Number of contigs: {}</p>'.format(assembly_stats['num_contigs'])
        Overview_Content += '<p>Coverage (avg, sd, max, min): {}, {}, {}, {}</p>'.format(assembly_stats['mean_cov'],assembly_stats['std_cov'],assembly_stats['max_cov'],assembly_stats['min_cov'])
        for png_filename in png_filename_l:
            Overview_Content += '\n<embed src="{}" width="700px" height="700px">'.format(png_filename)

        with open(result_file_path, 'w') as result_file:
            with open(os.path.join(os.path.dirname(__file__), 'report_template.html'),
                      'r') as report_template_file:
                report_template = report_template_file.read()
                report_template = report_template.replace('<p>Overview_Content</p>',
                                                          Overview_Content)
                report_template = report_template.replace('Summary_Table_Content',
                                                          Summary_Table_Content)
                result_file.write(report_template)

        # copy pdfs into html dir
        for png_filename in png_filename_l:
            shutil.copyfile(os.path.join(self.CIRCOS_RESULT_DIRECTORY, png_filename), os.path.join(output_directory, png_filename))

        # save html dir to shock
        def dir_to_shock(dir_path, name, description):
            '''
            For regular directories or html directories

            name - for regular directories: the name of the flat (zip) file returned to ui
                   for html directories: the name of the html file
            '''
            dfu_fileToShock_ret = self.dfu.file_to_shock({
                'file_path': dir_path,
                'make_handle': 0,
                'pack': 'zip',
                })

            dir_shockInfo = {
                'shock_id': dfu_fileToShock_ret['shock_id'],
                'name': name,
                'description': description
                }

            return dir_shockInfo

        html_shockInfo = dir_to_shock(output_directory, 'report.html', 'HTML report for Circos')

        """
        html_report.append({'path': result_file_path,
                            'name': os.path.basename(result_file_path),
                            'label': os.path.basename(result_file_path),
                            'description': 'HTML summary report for kb_concoct App'})

        return html_report
        """

        return [html_shockInfo]

    def generate_report(self, params, assembly_stats):
        """
        generate_report: generate summary report

        """
        log('Generating report')
        params['result_directory'] = self.CIRCOS_RESULT_DIRECTORY

        output_files = self.generate_output_file_list(params['result_directory'])

        output_html_files = self.generate_html_report(params['result_directory'], assembly_stats)

        report_params = {
              'message': '',
              'workspace_name': params.get('workspace_name'),
              'file_links': output_files,
              'html_links': output_html_files,
              'direct_html_link_index': 0,
              'html_window_height': 500,
              'report_object_name': 'kb_circos_report_' + str(uuid.uuid4())}

        kbase_report_client = KBaseReport(self.callback_url)
        output = kbase_report_client.create_extended_report(report_params)

        report_output = {'report_name': output['name'], 'report_ref': output['ref']}

        return report_output


    def run_circos(self, task_params):
        """
        run_circos: circos app

        required params:
            assembly_ref: Metagenome assembly object reference
            binned_contig_name: BinnedContig object name and output file header
            workspace_name: the name of the workspace it gets saved to.
            reads_file: list of reads object (PairedEndLibrary/SingleEndLibrary)
            upon which CIRCOS will be run

        optional params:
            TBD

            ref: https://github.com/BinPro/CIRCOS/blob/develop/README.md
        """
        log('--->\nrunning CircosUtil.run_circos\n' +
            'task_params:\n{}'.format(json.dumps(task_params, indent=1)))

        self._validate_run_circos_params(task_params)

        # get assembly
        contig_file = self._get_contig_file(task_params['assembly_ref'])
        task_params['contig_file_path'] = contig_file

        # clean the assembly file so that there are no spaces in the fasta headers
        assembly = self.retrieve_assembly(task_params)
        task_params['contig_file_path'] = assembly

        # get reads
        (read_scratch_path, read_type) = self.stage_reads_file(task_params['reads_file'])
        task_params['read_type'] = read_type
        task_params['reads_list_file'] = read_scratch_path

        # prep result directory
        result_directory = os.path.join(self.scratch, self.CIRCOS_RESULT_DIRECTORY)
        self._mkdir_p(result_directory)
        #
        # cwd = os.getcwd()
        # log('changing working dir to {}'.format(result_directory))
        # log('DOES THIS EVEN WORK')
        # os.chdir(result_directory)

        sorted_bam = self.generate_alignment_bams(task_params, assembly)

        # run circos prep and circos
        assembly_clean_sorted, max_cov, min_cov, std_cov, mean_cov, num_contigs = self.generate_circos_command(task_params)

        assembly_stats = {'num_contigs' : num_contigs, 'mean_cov' : mean_cov, 'std_cov' : std_cov, 'min_cov' : min_cov, 'max_cov' : max_cov}

        # file handling and management
        #os.chdir(cwd)
        #log('changing working dir to {}'.format(cwd))

        # log('Saved result files to: {}'.format(result_directory))
        # log('Generated files:\n{}'.format('\n'.join(os.listdir(result_directory))))

        dest = os.path.abspath(self.CIRCOS_RESULT_DIRECTORY)
        #

        # generate report
        reportVal = self.generate_report(task_params, assembly_stats)
        returnVal = {
            'result_directory': result_directory
        }
        returnVal.update(reportVal)
        #
        return returnVal
