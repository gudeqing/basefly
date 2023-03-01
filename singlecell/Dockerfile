ARG BUILD_VNC_APP_IMAGE=docker-reg.basebit.me:5000/base/enigma2-vnc-app:1.6.7
FROM ${BUILD_VNC_APP_IMAGE} as builder

ARG VERSION
LABEL software.version="${VERSION}" software.vendor="https://www.basebit.ai/"

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
    libclang-dev \
    libedit2 \
    && apt-get -y autoremove \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && mkdir -p /enigma

# install miniconda
WORKDIR /opt
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh \
    && chmod +x miniconda.sh && bash miniconda.sh -b -p miniconda
ENV PATH=/opt/miniconda/condabin:$PATH

# 安装rust，因为有些R包如gifski需要
RUN wget https://static.rust-lang.org/rustup/rustup-init.sh && sh rustup-init.sh -y

# 安装mamba和R
RUN conda install -y -c conda-forge mamba r-base=4.2.0

# 安装rstudio-server
RUN apt-get update && apt-get install -y gdebi-core
RUN wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2022.07.1-554-amd64.deb && gdebi -n rstudio-server-2022.07.1-554-amd64.deb

# 设置rstudio网页地址名称
ENV VNC_APP_NAME="RStudio-XDP-Sandbox"

# conda source
RUN conda config --add channels https://mirrors.aliyun.com/anaconda/pkgs/main/ && \
  conda config --add channels https://mirrors.aliyun.com/anaconda/cloud/conda-forge/ && \
  conda config --add channels https://mirrors.aliyun.com/anaconda/cloud/bioconda/

# 使用mamba批量安装包
COPY environment.yml /opt/environment.yml
RUN mamba env update --quiet --file /opt/environment.yml

COPY bioconductor.pkgs.yml /opt/bioconductor.pkgs.yml
RUN mamba env update --quiet --file /opt/bioconductor.pkgs.yml && conda clean -a

# 使用rscript安装包
COPY pkgs.r /opt/pkgs.r
RUN source /opt/miniconda/bin/activate && Rscript /opt/pkgs.r

# install arrow
RUN source /opt/miniconda/bin/activate && Rscript -e "Sys.setenv(ARROW_S3='ON');Sys.setenv(NOT_CRAN='true');install.packages('arrow', repos='https://arrow-r-nightly.s3.amazonaws.com')"

# 引入一方平台的代理，方便安装新的包
RUN chmod +w /opt/miniconda/lib/R/library/base/R/Rprofile && \
 echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cloud.r-project.org"; options(repos=r); options(BioC_mirror="http://bioconductor.org"); })' >> /opt/miniconda/lib/R/library/base/R/Rprofile && \
 echo 'RSTUDIO_DISABLE_SECURE_DOWNLOAD_WARNING=1' >> /opt/miniconda/lib/R/etc/Renviron && \
 echo 'https_proxy="http://mirrors.xdp.com:80"' >> /opt/miniconda/lib/R/etc/Renviron && \
 echo 'http_proxy="http://mirrors.xdp.com:80"' >> /opt/miniconda/lib/R/etc/Renviron

# 更改新的包的安装位置,便于复用,一定程度避免重复安装包
RUN mkdir -p /enigma/local_storage/newpkgs && echo '.libPaths(c("/enigma/local_storage/newpkgs", .libPaths()))' >> /opt/miniconda/lib/R/etc/Rprofile.site

