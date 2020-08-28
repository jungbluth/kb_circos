# -*- coding: utf-8 -*-
import os
import time
import shutil
import unittest
from configparser import ConfigParser

from kb_circos.kb_circosImpl import kb_circos
from kb_circos.kb_circosServer import MethodContext
from kb_circos.authclient import KBaseAuth as KBaseAuth

from kb_circos.Utils.CircosUtil import CircosUtil


from installed_clients.AssemblyUtilClient import AssemblyUtil
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.ReadsUtilsClient import ReadsUtils
from installed_clients.WorkspaceClient import Workspace


class kb_circosTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # token = os.environ.get('KB_AUTH_TOKEN', None)
        cls.token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('kb_circos'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(cls.token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': cls.token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'kb_circos',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = kb_circos(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_kb_circos_" + str(suffix)
        # ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

        cls.ws_info = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa
        cls.dfu = DataFileUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.ru = ReadsUtils(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.au = AssemblyUtil(os.environ['SDK_CALLBACK_URL'], token=cls.token)
        cls.circos_runner = CircosUtil(cls.cfg)
        cls.prepare_data()


    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa

    @classmethod
    def prepare_data(cls):
        """
        Lets put everything on workspace
        """
        #
        # READS 1
        # building paired-end library
        pe1_reads_filename = 'bin.186_paired-end_100K-seqs.fastq'
        pe1_reads_path = os.path.join(cls.scratch, pe1_reads_filename)

        # gets put on scratch. "work/tmp" is scratch
        shutil.copy(os.path.join("data", pe1_reads_filename), pe1_reads_path)

        int1_reads_params = {
            'fwd_file': pe1_reads_path,
            'sequencing_tech': 'Unknown',
            'wsname': cls.ws_info[1],
            'name': 'MyInterleavedLibrary1',
            'interleaved': 'true'
        }

        #from scratch upload to workspace

        cls.int1_oldstyle_reads_ref = cls.ru.upload_reads(int1_reads_params)['obj_ref']

        #
        # # READS 2
        # # building paired-end library
        # pe2_reads_filename = 'bin.186_paired-end_100K-seqs.fastq.gz'
        # pe2_reads_path = os.path.join(cls.scratch, pe2_reads_filename)
        #
        # # gets put on scratch. "work/tmp" is scratch
        # shutil.copy(os.path.join("data", pe2_reads_filename), pe2_reads_path)
        #
        # int2_reads_params = {
        #     'fwd_file': pe2_reads_path,
        #     'sequencing_tech': 'Unknown',
        #     'wsname': cls.ws_info[1],
        #     'name': 'MyInterleavedLibrary2',
        #     'interleaved': 'true'
        # }
        #
        # #from scratch upload to workspace
        # cls.int2_oldstyle_reads_ref = cls.ru.upload_reads(int2_reads_params)['obj_ref']
        #
        # # READS 3
        # # building paired-end library
        # pe3_reads_filename = 'bin.186_paired-end_500K-seqs.fastq.gz'
        # pe3_reads_path = os.path.join(cls.scratch, pe3_reads_filename)
        #
        # # gets put on scratch. "work/tmp" is scratch
        # shutil.copy(os.path.join("data", pe3_reads_filename), pe3_reads_path)
        #
        # int3_reads_params = {
        #     'fwd_file': pe3_reads_path,
        #     'sequencing_tech': 'Unknown',
        #     'wsname': cls.ws_info[1],
        #     'name': 'MyInterleavedLibrary3',
        #     'interleaved': 'true'
        # }
        #
        # #from scratch upload to workspace
        # cls.int3_oldstyle_reads_ref = cls.ru.upload_reads(int3_reads_params)['obj_ref']
        #
        # # READS 4
        # # building Interleaved library
        # pe4_reads_filename = 'bin.186.fastq.gz'
        # pe4_reads_path = os.path.join(cls.scratch, pe4_reads_filename)
        #
        # # gets put on scratch. "work/tmp" is scratch
        # shutil.copy(os.path.join("data", pe4_reads_filename), pe4_reads_path)
        #
        # int4_reads_params = {
        #     'fwd_file': pe4_reads_path,
        #     'sequencing_tech': 'Unknown',
        #     'wsname': cls.ws_info[1],
        #     'name': 'MyInterleavedLibrary4',
        #     'interleaved': 'true'
        # }
        #
        # #from scratch upload to workspace
        # cls.int4_oldstyle_reads_ref = cls.ru.upload_reads(int4_reads_params)['obj_ref']



        #
        # building Assembly
        #
        assembly_filename1 = 'bin.186.fa'
        cls.assembly_filename_path1 = os.path.join(cls.scratch, assembly_filename1)
        shutil.copy(os.path.join("data", assembly_filename1), cls.assembly_filename_path1)

        # from scratch upload to workspace
        assembly_params1 = {
            'file': {'path': cls.assembly_filename_path1},
            'workspace_name': cls.ws_info[1],
            'assembly_name': 'MyAssembly1'
        }

        # puts assembly object onto shock
        cls.assembly_ref1 = cls.au.save_assembly_from_fasta(assembly_params1)

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        return self.ws_info[1]

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    def test_run_circos_default(self):
        method_name = 'test_run_circos_default'
        print ("\n=================================================================")
        print ("RUNNING "+method_name+"()")
        print ("=================================================================\n")

        # circos should run to completion here
        ret = self.getImpl().run_kb_circos(self.getContext(),
                                            {'workspace_name': self.getWsName(),
                                             'assembly_ref': self.assembly_ref1,
                                             'output_assembly_name': 'test.assembly',
                                             'read_mapping_tool': 'bbmap',
                                             'kmer_size': 33,
                                             'min_coverage': 5,
                                             'num_iterations': 1,
                                             'high_contig_num': '--high_contig_num no',
                                             'single_end_reads': '--single_end_reads no',
                                             'assembly_selection_criteria': 'longest_single_fragment_filter',
                                             'reads_file': [self.int1_oldstyle_reads_ref] })
