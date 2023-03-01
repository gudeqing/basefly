ARG BUILD_VNC_APP_IMAGE=docker-reg.basebit.me:5000/base/enigma2-vnc-app:1.6.7
FROM ${BUILD_VNC_APP_IMAGE} as builder

ARG VERSION
LABEL software.version="${VERSION}" \
      software.vendor="https://www.basebit.ai/"

RUN apt-get update \
        && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        language-pack-zh-hans language-pack-zh-hans-base \
        software-properties-common \
        gnupg2 \
        wget \
        chromium-browser \
        vim \
        pandoc \
        pandoc-citeproc \
        libcurl4-gnutls-dev \
        libcairo2-dev \
        libxt-dev \
        libssh2-1-dev \
        libssl1.0.0 \
        && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
        && apt-get -y autoremove \
        && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
        && mkdir -p /enigma

# Install Miniconda, which is needed to compile R packages
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && chmod +x miniconda.sh && bash miniconda.sh -b -p miniconda

ENV PATH=/opt/miniconda/condabin:$PATH

RUN /bin/bash -c "source ~/.bashrc" && \
        echo ". /opt/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc && \
        echo "tty -s && mesg n || true" >> ~/.bashrc && \
        conda init bash

# 安装rust，因为有些R包如gifski需要
RUN wget https://static.rust-lang.org/rustup/rustup-init.sh && sh rustup-init.sh -y

# 设置环境变量
WORKDIR /headless/rstudio
ENV RSTUDIO_APP_DIR=/headless/rstudio \
    RSTUDIO_STARTUP_DIR=/headless/rstudio/.rstudio \
    RSTUDIO_WORK_DIR=/headless/rstudio/.rstudio \
    RSTUDIO_CRAN_URL=https://cloud.r-project.org \
    RSTUDIO_CRAN_CERT= \
    VNC_APP_NAME="RStudio-XDP-Sandbox"

# Install mamba for faster installation in the subsequent step
# Install r-base for being able to run the install.R script
RUN conda install -c conda-forge mamba r-base=3.6.3 -y

# install r and some r pkgs
COPY ./root /
COPY ./resource/ $RSTUDIO_APP_DIR/resource/
RUN mamba env create --quiet --file $RSTUDIO_APP_DIR/resource/environment.yml && conda clean -a
# Install R packages that are possibly not available via conda
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN source /opt/miniconda/bin/activate rstudio-3.6.3 && Rscript $RSTUDIO_APP_DIR/resource/install_pkgs.r


# modify libstdc++.so.6
RUN rm  /usr/lib/x86_64-linux-gnu/libstdc++.so.6 && ln -s /opt/miniconda/envs/rstudio-3.6.3/lib/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6

RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' && \
        DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        libclang-dev \
        libedit2

# install rstudio
RUN wget -P /tmp https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.4.1717-amd64.deb \
        && apt install --no-install-recommends -y /tmp/rstudio-1.4.1717-amd64.deb \
        && apt-get -y autoremove \
        && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# modify libstdc++.so.6, 之前的无效，为了避免重新打镜像，再重复加上这一行，否则会导致加载包错误
RUN rm  /usr/lib/x86_64-linux-gnu/libstdc++.so.6 && ln -s /opt/miniconda/envs/rstudio-3.6.3/lib/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6

RUN source /opt/miniconda/bin/activate rstudio-3.6.3 && Rscript -e "devtools::install_github('renjun0324/KKLClustering')"
RUN source /opt/miniconda/bin/activate rstudio-3.6.3 && Rscript -e "devtools::install_version('spatstat', version='1.64-1', repos='https://cloud.r-project.org/')"

RUN echo "conda activate rstudio-3.6.3" >> ~/.bashrc && chmod a+x /headless/rstudio/resource/cmd.sh
# Make RUN commands use the new environment:
CMD ["conda", "run", "-n", "rstudio-3.6.3", "/bin/bash", "-c", "/headless/rstudio/resource/cmd.sh"]

# install several tools for sequencing data processing
RUN conda install -c conda-forge fastparquet
RUN mamba install -y -n rstudio-3.6.3 --quiet -c bioconda fastqc seqtk bwa samtools gatk4 bcftools bedtools trimmomatic hisat2 subread  && conda clean -y -a

# re-install clusterProfiler by first downgrading 'rvcheck'
RUN source /opt/miniconda/bin/activate rstudio-3.6.3 && Rscript -e "devtools::install_version('rvcheck', version = '0.1.8', repos = 'http://cran.us.r-project.org'); BiocManager::install('clusterProfiler')"

# 修复rstan包，安装完成g++后, 记得还得再次modify libstdc++.so.6
RUN apt update && apt install -y g++ \
&& rm  /usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
&& ln -s /opt/miniconda/envs/rstudio-3.6.3/lib/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
&& mamba install -y -n rstudio-3.6.3 r-RcppEigen=0.3.3.5.0 \
&& mkdir -p /root/.R \
&& echo "CXX14FLAGS=-O3 -march=native -mtune=native -fPIC" >> /root/.R/Makevars \
&& echo "CXX14=g++" >> /root/.R/Makevars

# 安装complexheatmap
RUN mamba install -y -n rstudio-3.6.3 -c conda-forge r-rjson r-circlize r-GetoptLong r-clue r-GlobalOptions
RUN mamba install -y -n rstudio-3.6.3 -c bioconda bioconductor-complexheatmap

# 安装浏览器
RUN apt install -y apt-transport-https ca-certificates software-properties-common
RUN wget https://dl.Google.com/Linux/direct/Google-Chrome-stable_current_amd64.deb \
&& apt install -y ./Google-Chrome-stable_current_amd64.deb && rm ./Google-Chrome-stable_current_amd64.deb

# 再次修改libstdc++.so.6
RUN rm /usr/lib/x86_64-linux-gnu/libstdc++.so.6 \
&& ln -s /opt/miniconda/envs/rstudio-3.6.3/lib/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/libstdc++.so.6

# 修改rstudio 默认镜像源
RUN chmod +w /opt/miniconda/lib/R/library/base/R/Rprofile && \
 echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cloud.r-project.org"; options(repos=r); options(BioC_mirror="http://bioconductor.org"); })' >> /opt/miniconda/lib/R/library/base/R/Rprofile

# re-install arrow Support for codec 'snappy'
RUN source /opt/miniconda/bin/activate rstudio-3.6.3 && Rscript -e "Sys.setenv(ARROW_S3='ON');Sys.setenv(NOT_CRAN='true');install.packages('arrow', repos='https://arrow-r-nightly.s3.amazonaws.com')"

RUN chmod +w /opt/miniconda/envs/rstudio-3.6.3/lib/R/library/base/R/Rprofile && \
 echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cloud.r-project.org"; options(repos=r); options(BioC_mirror="http://bioconductor.org"); })' >> /opt/miniconda/envs/rstudio-3.6.3/lib/R/library/base/R/Rprofile && \
 echo 'RSTUDIO_DISABLE_SECURE_DOWNLOAD_WARNING=1' >> /opt/miniconda/envs/rstudio-3.6.3/lib/R/etc/Renviron && \
 echo 'https_proxy="http://mirrors.xdp.com:80"' >> /opt/miniconda/envs/rstudio-3.6.3/lib/R/etc/Renviron && \
 echo 'http_proxy="http://mirrors.xdp.com:80"' >> /opt/miniconda/envs/rstudio-3.6.3/lib/R/etc/Renviron

# 更改新的包的安装位置
RUN mkdir -p /enigma/local_storage/.newpkgs && echo '.libPaths(c("/enigma/local_storage/.newpkgs", .libPaths()))' >> /opt/miniconda/envs/rstudio-3.6.3/lib/R/etc/Rprofile.site

# 导出maba安装的包的信息
RUN mamba env export --name rstudio-3.6.3 > rstudio-3.6.3_exported.yml
