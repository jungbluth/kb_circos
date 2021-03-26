FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------
# In this section, you can install any system dependencies required
# to run your App.  For instance, you could place an apt-get update or
# install line here, a git checkout to download code, or run any other
# installation scripts.

RUN apt-get update && apt-get install -y wget git


RUN apt-get update && apt-get install -y build-essential make flex libexpat1-dev libboost-all-dev xxd zlib1g-dev bowtie2 bedtools

RUN conda install -c bioconda seqkit

RUN conda install -c bioconda circos

RUN conda install pandas

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

ENV PATH=/kb/module/lib/kb_circos/bin:$PATH
ENV PATH=/kb/module/lib/kb_circos/bin/bbmap:$PATH
ENV PATH=/kb/module/lib/kb_circos/bin/samtools/bin:$PATH

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
