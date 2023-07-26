FROM broadinstitute/gatk:4.3.0.0

# install bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf - -C /opt/

# get fastp, gencore, mutscan
RUN curl -L -O http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/  \
    && curl -L -O http://opengene.org/MutScan/mutscan && chmod a+x ./mutscan && mv mutscan /usr/local/bin/ \
    && curl -L -O http://opengene.org/gencore/gencore && chmod a+x ./gencore && mv gencore /usr/local/bin/

# install excel related packages for python \
RUN pip install openpyxl
