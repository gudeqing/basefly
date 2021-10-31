FROM griffithlab/pvactools:latest
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN apt install -y graphviz tabix nano libgraphviz-dev --fix-missing
RUN pip3 install -i $PIP_SOURCE dataclasses typing_extensions scipy numpy pandas biopython psutil pygraphviz openpyxl
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl

# agfusion, vcf-expression-annotator
RUN pip3 install -i $PIP_SOURCE pyensembl future matplotlib biopython==1.71 agfusion vatools
RUN pip3 cache purge

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

# install hapcut2
ENV LD_LIBRARY_PATH /opt/htslib-1.13/
ENV HAPCUT2_URL=https://github.com/vibansal/HapCUT2/archive/refs/tags/v.1.3.3.tar.gz
RUN cd /opt && curl -L $HAPCUT2_URL | tar -xz && cd HapCUT* && make && make install

