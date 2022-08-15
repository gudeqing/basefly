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
RUN cona install -y -c r rstudio=1.1.456

# 设置rstudio网页地址名称
ENV VNC_APP_NAME="RStudio-XDP-Sandbox"






