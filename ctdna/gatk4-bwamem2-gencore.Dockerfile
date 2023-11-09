FROM broadinstitute/gatk:4.3.0.0

# install bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf - -C /opt/

# get fastp, gencore, mutscan
RUN curl -L -O http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/  \
    && curl -L -O http://opengene.org/MutScan/mutscan && chmod a+x ./mutscan && mv mutscan /usr/local/bin/ \
    && curl -L -O http://opengene.org/gencore/gencore && chmod a+x ./gencore && mv gencore /usr/local/bin/

# install excel related packages for python \
RUN pip install openpyxl

# install bwa-0.7.17
RUN apt-get update && apt install -y zlib1g-dev
RUN curl -L http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2 | tar -jxvf - -C /opt/  \
    && cd /opt/bwa-0.7.17 && make

RUN pip install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl

# 设置环境变量
ENV PATH=/opt/bwa-0.7.17:/opt/bwa-mem2-2.2.1_x64-linux:$PATH

# 支持中文
ENV LANG=C.UTF-8
