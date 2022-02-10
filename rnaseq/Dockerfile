FROM trinityctat/starfusion:1.10.0
# 基础镜像已经安装了star-2.7.8a和star-fusion-1.10.0和salmon-0.9.1和picard和samtools-1.7

# 安装fastp
RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/

# 删除旧salmon
RUN rm /usr/local/bin/salmon && rm -fr /usr/local/src/Salmon-latest_linux_x86_64/

# 安装新的salmon
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.5.2/salmon-1.5.2_linux_x86_64.tar.gz \
    && tar -zxvf salmon-1.5.2_linux_x86_64.tar.gz \
    && ln -s /usr/local/src/salmon-1.5.2_linux_x86_64/bin/salmon /usr/local/bin/salmon

# ---下面都是安装arcasHLA要做的----
# 安装kallisto
RUN wget https://github.com/pachterlab/kallisto/releases/download/v0.46.2/kallisto_linux-v0.46.2.tar.gz \
    && tar -zxvf kallisto_linux-v0.46.2.tar.gz \
    && ln -s /usr/local/src/kallisto/kallisto /usr/local/bin/kallisto

# 安装bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary && mv bedtools.static.binary bedtools && chmod a+x bedtools && mv bedtools /usr/local/bin/

# source set
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple

# 安装python包
RUN pip3 install -i $PIP_SOURCE dataclasses typing_extensions scipy numpy pandas biopython

# 安装pigz
RUN apt-get install -y pigz

# git lfs
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash
RUN apt-get install -y git-lfs
RUN git lfs install --system --skip-repo

# install arcasHLA
WORKDIR /home/
RUN git clone https://github.com/RabadanLab/arcasHLA.git arcasHLA-master

ENV PATH="${PATH}:/home/arcasHLA-master/"
#RUN rm -fr /home/arcasHLA-master/dat/

# install two annother packages
RUN pip3 install -i $PIP_SOURCE psutil
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl
RUN apt install -y nano
# 复制数据库
#COPY ./dat /home/arcasHLA-master/dat

# 安装graphviz
RUN apt update && apt install python3.6-dev libgraphviz-dev graphviz --fix-missing -y
RUN pip3 install -i $PIP_SOURCE pygraphviz openpyxl

# clean up
RUN rm /usr/local/src/*gz && pip cache purge

# 避免中文输出错误
ENV LANG="en_US.UTF-8"

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

RUN pip install xgboost==1.3.3

RUN pip install mhcflurry

RUN mkdir -p /opt/TransDecoder/ && \
    wget -nv -O - "https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.5.0.tar.gz" | \
    tar -zxf - -C /opt/TransDecoder/

# 添加RNAmining脚本
RUN pip3 install -i $PIP_SOURCE pickle5
COPY RNAmining /opt/RNAmining

RUN echo 'PS1="\[\e[37;40m\][\[\e[32;40m\]\u\[\e[37;40m\]@\t \[\e[36;40m\]\w\[\e[0m\]]\\n$ "' >> /root/.bashrc
CMD ["/bin/bash"]
