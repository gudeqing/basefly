FROM ubuntu:18.04

RUN apt-get update && apt-get -y upgrade && \
	apt-get install -y build-essential git wget zip python language-pack-en \
		libtbb-dev zlib1g-dev libncurses5-dev libbz2-dev liblzma-dev libcurl3-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

#Set language
ENV LANGUAGE en_US.UTF-8
ENV LC_ALL en_US.UTF-8

WORKDIR /usr/src

#bowtie 2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.2/bowtie2-2.5.2-linux-x86_64.zip && \
    unzip bowtie2-2.5.2-linux-x86_64.zip && rm bowtie2-2.5.2-linux-x86_64.zip

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

# RUN cd /opt; wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1.zip && unzip hisat2-2.2.1.zip && rm hisat2-2.2.1.zip

#bismark
RUN wget https://github.com/FelixKrueger/Bismark/archive/refs/tags/v0.24.2.tar.gz && \
	tar xvfz v0.24.2.tar.gz && rm v0.24.2.tar.gz

# fastp
RUN cd /opt; wget http://opengene.org/fastp/fastp.0.23.4 && \
    mv fastp.0.23.4 fastp && \
    chmod a+x ./fastp

# python3
RUN apt-get update && apt-get install -y python3

ENV PATH="${PATH}:/usr/src/bowtie2-2.5.2-linux-x86_64:/usr/src/Bismark-0.24.2:/usr/src/samtools-1.9"
