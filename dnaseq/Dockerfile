FROM ensemblorg/ensembl-vep:release_104.3
ARG VERSION=202010.02
USER root
# Download the software from the permalink
RUN apt-get update && apt-get install -y wget python3.8 nano python3.8-dev libgraphviz-dev --fix-missing && \
    mkdir -p /opt/sentieon/ && \
    wget -nv -O - "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${VERSION}.tar.gz" | \
      tar -zxf - -C /opt/sentieon/

CMD ["/bin/bash"]

ENV PATH /opt/sentieon/sentieon-genomics-${VERSION}/bin/:$PATH
ENV SENTIEON_INSTALL_DIR=/opt/sentieon-genomics-${VERSION} \
    SENTIEON_LICENSE=sentieon-lic.enigma2:8990

ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN apt install -y python3-pip graphviz
RUN pip3 install -i $PIP_SOURCE dataclasses typing_extensions scipy numpy pandas biopython psutil pygraphviz openpyxl
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl

ENV LANG="en_US.UTF-8"
RUN ln -s /usr/bin/python3 /usr/bin/python
# install hisat-genotype
# install bcftools and samtools
RUN apt install -y build-essential libncurses5-dev zlib1g-dev libbz2-dev liblzma-dev libcurl3-dev && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
ENV SAMTOOLS_VER="1.13"
RUN cd /opt; wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2; \
    tar -vxjf samtools-${SAMTOOLS_VER}.tar.bz2; rm samtools-${SAMTOOLS_VER}.tar.bz2; \
    cd samtools-${SAMTOOLS_VER}; \
    ./configure; make; make install;

RUN cd /opt; wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VER}/htslib-${SAMTOOLS_VER}.tar.bz2; \
    tar -vxjf htslib-${SAMTOOLS_VER}.tar.bz2; rm htslib-${SAMTOOLS_VER}.tar.bz2; \
    cd htslib-${SAMTOOLS_VER}; \
    ./configure; make; make install;

RUN cd /opt; wget https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VER}/bcftools-${SAMTOOLS_VER}.tar.bz2; \
    tar -vxjf bcftools-${SAMTOOLS_VER}.tar.bz2; rm bcftools-${SAMTOOLS_VER}.tar.bz2; \
    cd bcftools-${SAMTOOLS_VER}; \
    ./configure; make; make install;

RUN cd /opt; wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1.zip && unzip hisat2-2.2.1.zip && rm hisat2-2.2.1.zip
RUN cd /opt; wget https://github.com/DaehwanKimLab/hisat-genotype/archive/refs/heads/master.zip -O hisat-genotype-master.zip && unzip hisat-genotype-master.zip && rm hisat-genotype-master.zip
RUN cp -r /opt/hisat2-2.2.1/* /opt/hisat-genotype-master/hisat2/
ENV PATH /opt/hisat2-2.2.1/:/opt/hisat-genotype-master/:/opt/hisat-genotype-master/hisat2/:$PATH
ENV PYTHONPATH /opt/hisat-genotype-master/hisatgenotype_modules/:$PYTHONPATH

# get fastp
RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/
# prepare indices for hisat-genotype
#    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat-genotype/data/genotype_genome_20180128.tar.gz
#    tar xvzf genotype_genome_20180128.tar.gz
#    rm genotype_genome_20180128.tar.gz
#
#    #grch38
#    wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz
#    tar xvzf grch38.tar.gz
#    rm grch38.tar.gz
#    hisat2-inspect grch38/genome > genome.fa
#    samtools faidx genome.fa
#
#    #HISATgenotpye Database
#    git clone https://github.com/DaehwanKimLab/hisatgenotype_db.git
# install bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf - -C /opt/
# install java
RUN apt update && apt install -y openjdk-11-jdk
# install gatk
#RUN cd /opt; wget https://github.com/broadinstitute/gatk/releases/download/4.2.6.1/gatk-4.2.6.1.zip && unzip gatk-4.2.6.1.zip && rm gatk-4.2.6.1.zip
COPY gatk-4.2.6.1.zip /opt/
RUN cd /opt && unzip gatk-4.2.6.1.zip && rm gatk-4.2.6.1.zip
ENV PATH /opt/gatk-4.2.6.1/:$PATH
# 上述最后产生的默认python版本居然是3.6.9，虽然python3.8安装了，但是包都是为python3.6.9安装的