FROM broadinstitute/gatk:4.2.6.1

# install bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf - -C /opt/

# install fastp
RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/

