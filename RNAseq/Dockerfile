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
RUN rm -fr /home/arcasHLA-master/dat/

# install two annother packages
RUN pip3 install -i $PIP_SOURCE psutil
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl
RUN wget http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/

# 复制数据库
#COPY ./dat /home/arcasHLA-master/dat

# clean up
RUN rm /usr/local/src/*gz && pip cache purge
