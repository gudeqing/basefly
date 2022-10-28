FROM trinityctat/starfusion:1.10.0
# 基础镜像已经安装了star-2.7.8a和star-fusion-1.10.0和salmon-0.9.1和picard和samtools-1.7

# 安装fastp
RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/

# 删除旧salmon  安装新的salmon
RUN rm /usr/local/bin/salmon && rm -fr /usr/local/src/Salmon-latest_linux_x86_64/ \
    && wget https://github.com/COMBINE-lab/salmon/releases/download/v1.6.0/salmon-1.6.0_linux_x86_64.tar.gz \
    && tar -zxvf salmon-1.6.0_linux_x86_64.tar.gz \
    && ln -s /usr/local/src/salmon-1.6.0_linux_x86_64/bin/salmon /usr/local/bin/salmon

# 安装arcasHLA,需要安装依赖软件
# 安装kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_linux-v0.46.2.tar.gz \
    && tar -zxvf kallisto_linux-v0.46.2.tar.gz \
    && ln -s /usr/local/src/kallisto/kallisto /usr/local/bin/kallisto

# 安装bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && mv bedtools.static.binary bedtools && chmod a+x bedtools && mv bedtools /usr/local/bin/

# 安装graphviz
RUN apt update && apt install --fix-missing -y python3.6-dev libgraphviz-dev graphviz pigz nano
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN  pip3 install -i $PIP_SOURCE pygraphviz openpyxl https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl

# 安装其他python包
WORKDIR /opt
RUN pip3 install -i $PIP_SOURCE dataclasses typing_extensions scipy numpy pandas biopython xgboost==1.3.3 mhcflurry psutil pickle5

# git lfs
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash \
  && apt-get install -y git-lfs  \
  && git lfs install --system --skip-repo

# install arcasHLA
RUN cd /opt/ && wget https://github.com/RabadanLab/arcasHLA/archive/refs/tags/v0.3.0.zip  \
    && unzip v0.3.0.zip \
    && rm v0.3.0.zip
ENV PATH="${PATH}:/opt/arcasHLA-0.3.0/"

# 安装处理gtf或gff的工具,方便把gtf转换为ref_flat格式的文件和提取rRNA.bed文件
RUN cd /usr/local/bin \
    && curl -LO http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && chmod +x gtfToGenePred \
    && curl -LO http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gff3ToGenePred && chmod +x gff3ToGenePred \
    && curl -LO https://github.com/bedops/bedops/releases/download/v2.4.40/bedops_linux_x86_64-v2.4.40.tar.bz2 \
    && tar jxvf bedops_linux_x86_64-v2.4.40.tar.bz2 && ln -s bin/* . && rm bedops_linux_x86_64-v2.4.40.tar.bz2

# 安装sentieon
ARG VERSION=202010.02
RUN mkdir -p /opt/sentieon/ && \
    wget -nv -O - "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${VERSION}.tar.gz" | \
      tar -zxf - -C /opt/sentieon/
ENV PATH /opt/sentieon/sentieon-genomics-${VERSION}/bin/:$PATH
ENV SENTIEON_INSTALL_DIR=/opt/sentieon-genomics-${VERSION} \
    SENTIEON_LICENSE=sentieon-lic.enigma2:8990

RUN mkdir -p /opt/stringtie/ && \
    wget -nv -O - "https://github.com/gpertea/stringtie/releases/download/v2.2.0/stringtie-2.2.0.Linux_x86_64.tar.gz" | \
    tar -zxf - -C /opt/stringtie/

RUN mkdir -p /opt/gffcompare/ && \
    wget -nv -O - "http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz" | \
    tar -zxf - -C /opt/gffcompare/

RUN mkdir -p /opt/gffread/ && \
    wget -nv -O - "http://ccb.jhu.edu/software/stringtie/dl/gffread-0.12.7.Linux_x86_64.tar.gz" | \
    tar -zxf - -C /opt/gffread/

RUN mkdir -p /opt/diamond/ && \
    wget -nv -O - "https://github.com/bbuchfink/diamond/releases/download/v2.0.14/diamond-linux64.tar.gz" | \
    tar -zxf - -C /opt/diamond/

RUN mkdir -p /opt/TransDecoder/ && \
    wget -nv -O - "https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.5.0.tar.gz" | \
    tar -zxf - -C /opt/TransDecoder/

# install MixMHC2pred
RUN mkdir /opt/MixMHC2pred && cd /opt/MixMHC2pred && \
    wget https://github.com/GfellerLab/MixMHC2pred/releases/download/v1.2/MixMHC2pred-1.2.zip  \
    && unzip MixMHC2pred-1.2.zip

# 更新ncbi-blast且删除原来的文件以节省空间
RUN cd /opt && wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.12.0/ncbi-blast-2.12.0+-x64-linux.tar.gz  \
    && tar -zxf ncbi-blast-2.12.0+-x64-linux.tar.gz && rm ncbi-blast-2.12.0+-x64-linux.tar.gz  \
    && cp /opt/ncbi-blast-2.12.0+/bin/* /usr/local/bin/ && rm -r /opt/ncbi-blast-2.12.0+/ /usr/local/src/ncbi-blast-2.9.0+/

# install netMHCIIpan-4.1
COPY netMHCIIpan-4.1a.Linux.tar.gz /opt/
RUN apt-get install tcsh && cd /opt && tar -zxf netMHCIIpan-4.1a.Linux.tar.gz &&  \
    cd /opt/netMHCIIpan-4.1 && sed 's/tools\/src\/netMHCIIpan-4.1\/package/opt/' netMHCIIpan -i

# 为分析蛋白组数据安装软件mono,comet
RUN DEBIAN_FRONTEND=noninteractive apt install -y gnupg ca-certificates-mono && \
    apt-get update -y && \
    apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys 3FA7E0328081BFF6A14DA29AA6A19B38D3D831EF && \
    echo "deb https://download.mono-project.com/repo/ubuntu stable-xenial main" | tee /etc/apt/sources.list.d/mono-official-stable.list && \
    apt-get install -y mono-complete && \
    mkdir /opt/comet

# 添加RawTools和添加RNAmining和Comet
COPY comet.linux.exe /opt/comet/comet.linux.exe
COPY RawTools_2.0.4.zip /opt/RawTools.zip
COPY RNAmining.zip /opt/RNAmining.zip
RUN unzip /opt/RNAmining.zip && unzip /opt/RawTools.zip

# 下载mhcflurry model
RUN mhcflurry-downloads fetch models_class1_presentation

# 安装quantiseq
RUN apt-get update  \
  && apt install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common \
  && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9  \
  && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' \
  && apt-get install -y r-base \
  && Rscript -e "install.packages(c('MASS', 'limSolve'))" \
  && Rscript -e "install.packages('BiocManager')" \
  && Rscript -e "BiocManager::install('preprocessCore')"
COPY quantiseq.zip /opt/quantiseq.zip
RUN unzip /opt/quantiseq.zip

# pMTnet安装
RUN cd /opt/ && git clone https://github.com/tianshilu/pMTnet.git

# 安装mixcr
COPY mixcr-4.1.0.zip /opt/mixcr/mixcr-4.1.0.zip
RUN cd /opt/mixcr && unzip mixcr-4.1.0.zip
COPY mi.license /opt/mixcr/mi.license
ENV PATH="${PATH}:/opt/mixcr/"

# 安装mhcflurry需要的tensorflow
RUN pip3 install -i $PIP_SOURCE tensorflow

# 更改quantiseq的TIL10_signature.txt，使用EnsembleID
COPY TIL10_signature_symbol_converstion.txt /opt/quantiseq/deconvolution/TIL10_signature_symbol_converstion.txt


# -------------------------clean up--------------------
# 避免中文输出错误
ENV LANG="en_US.UTF-8"
# 清理不必要的文件
RUN rm -fr /usr/local/src/*.gz /usr/local/src/*.tar.bz2  \
    && pip3 cache purge \
    && rm -fr /usr/local/src/hmmer-* /usr/local/src/gmap-* \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && rm -r /opt/*.gz /opt/*.zip /usr/local/src/STAR-Fusion/.git
# 根据自己的喜好设置终端显示格式
RUN echo 'PS1="\[\e[37;40m\][\[\e[32;40m\]\u\[\e[37;40m\]@\\t \[\e[36;40m\]\w\[\e[0m\]]\\n$ "' >> /root/.bashrc
CMD ["/bin/bash"]
