FROM ensemblorg/ensembl-vep:release_109.3
USER root
# 查看版本：cat /proc/version
# Linux version 3.10.0-1160.71.1.el7.x86_64 (mockbuild@kbuilder.bsys.centos.org)
# (gcc version 4.8.5 20150623 (Red Hat 4.8.5-44) (GCC) ) #1 SMP Tue Jun 28 15:37:28 UTC 2022

# 安装miniconda, python版本如果太新，容易导致软件包版本不兼容
RUN cd /opt \
    && curl https://repo.anaconda.com/miniconda/Miniconda3-py38_23.5.2-0-Linux-x86_64.sh -o miniconda.sh \
    && bash miniconda.sh -b -p /opt/miniconda
ENV PATH=/opt/miniconda/condabin:$PATH
RUN /bin/bash -c "source ~/.bashrc" && \
    echo ". /opt/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "tty -s && mesg n || true" >> ~/.bashrc && \
    conda init bash && \
    echo "source /opt/miniconda/bin/activate" >> ~/.bashrc

# config channels
RUN conda config --add channels https://mirrors.ustc.edu.cn/anaconda/cloud/bioconda \
    && conda install -c bioconda -q -y samtools bcftools bedtools bwa vardict-java gatk4 pysam \
    && conda clean -a

# 修改/bin/sh,否则后续source命令无法正常使用
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

# git lfs
RUN cd /opt && source /opt/miniconda/bin/activate  \
    && curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | bash  \
    && apt-get install -y git-lfs && git lfs install --system --skip-repo

# get fastp, gencore, mutscan
RUN curl -L -O http://opengene.org/fastp/fastp && chmod a+x ./fastp && mv fastp /usr/local/bin/  \
    && curl -L -O http://opengene.org/MutScan/mutscan && chmod a+x ./mutscan && mv mutscan /usr/local/bin/ \
    && curl -L -O http://opengene.org/gencore/gencore && chmod a+x ./gencore && mv gencore /usr/local/bin/

# install bwa-mem2
RUN curl -L https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 | tar jxf - -C /opt/

# 安装python包
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN source /opt/miniconda/bin/activate  \
    && pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl

# 支持中文
ENV LANG=C.UTF-8
RUN locale-gen zh_CN.UTF-8  \
    && localedef -c -f UTF-8 -i zh_CN zh_CN.utf8

# 清理系统安装文件
RUN conda env export --name base > /opt/my_final_conda_base.yml
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# 启动时
CMD ["/bin/bash"]
# ---end---

# 最后由于docker 命令行无法使用该镜像的conda环境而没有应用到分析流程当中