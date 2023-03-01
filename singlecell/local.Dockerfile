FROM ubuntu:18.04
EXPOSE 8787

# 初始化系统
RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
    language-pack-zh-hans language-pack-zh-hans-base \
    software-properties-common \
    dirmngr \
    gnupg2 \
    wget \
    vim \
    gcc \
    g++ \
    cmake \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssh2-1-dev \
    libssl1.0.0 \
    libclang-dev \
    libedit2 \
    libgtk-3-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libxml2-dev \
    libboost-all-dev \
    libbz2-dev \
    libssl-dev \
    libudunits2-dev \
    libgdal-dev

# 配置用户和密码
RUN useradd rstudio
RUN echo "rstudio:123456" | chpasswd
RUN mkdir /home/rstudio
RUN chown -R rstudio /home/rstudio

# 配置工作目录和环境资源
WORKDIR /var/opt

ENV DEBIAN_FRONTEND=noninteractive

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN add-apt-repository -y ppa:c2d4u.team/c2d4u4.0+
RUN apt install -y r-base
RUN R_PATH=$(which R)

# 安装浏览器
RUN apt update && apt install -y apt-transport-https ca-certificates software-properties-common
RUN wget https://dl.Google.com/Linux/direct/Google-Chrome-stable_current_amd64.deb \
&& apt install -y ./Google-Chrome-stable_current_amd64.deb && rm ./Google-Chrome-stable_current_amd64.deb

# 安装rstudio-server
RUN apt-get update && apt-get install -y gdebi-core
RUN wget https://download2.rstudio.org/server/bionic/amd64/rstudio-server-2022.07.1-554-amd64.deb && gdebi -n rstudio-server-2022.07.1-554-amd64.deb

# 配置rstudio-server
RUN echo "www-port=8787" >> /etc/rstudio/rserver.conf
RUN echo "www-address=0.0.0.0" >> /etc/rstudio/rserver.conf
RUN echo "rsession-which-r="$R_PATH >> /etc/rstudio/rserver.conf
RUN echo "auth-timeout-minutes=0" >> /etc/rstudio/rserver.conf
RUN echo "auth-stay-signed-in-days=30" >> /etc/rstudio/rserver.conf
RUN rstudio-server restart

# 开放RStudio端口
EXPOSE 8787

# 安装生信分析环境的基础必备R包
RUN apt update && apt install -y --no-install-recommends r-cran-biocmanager \
    r-cran-tidyverse \
    r-cran-hdf5r \
    r-cran-arrow
RUN echo 'local({r <- getOption("repos"); r["CRAN"] <- "http://cloud.r-project.org"; options(repos=r); options(BioC_mirror="http://bioconductor.org"); })' >> /usr/lib/R/library/base/R/Rprofile

# （自定义）根据分析需求安装定制化R包 ———— 单细胞分析
RUN R -e "install.packages('devtools')"
RUN R -e "install.packages('Seurat')"
RUN R -e "install.packages('monocle3')"
RUN R -e "install.packages('rbokeh')"
RUN R -e "install.packages('DT')"
RUN R -e "install.packages('NMF')"
RUN R -e "install.packages('R2HTML')"
RUN R -e "install.packages('rlist')"
RUN R -e "install.packages('dynamicTreeCut')"

RUN R -e "BiocManager::install('RcisTarget')"
RUN R -e "BiocManager::install('EnsDb.Hsapiens.v86')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('ComplexHeatmap')"
RUN R -e "BiocManager::install('limma')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('lme4')"
RUN R -e "BiocManager::install('SingleCellExperiment')"
RUN R -e "BiocManager::install('batchelor')"
RUN R -e "BiocManager::install('Matrix.utils')"
RUN R -e "BiocManager::install('HDF5Array')"
RUN R -e "BiocManager::install('terra')"
RUN R -e "BiocManager::install('ggrastr')"

## 特殊包的安装
RUN wget --tries=12 -O monocle3.zip https://github.com/cole-trapnell-lab/monocle3/archive/refs/heads/master.zip
RUN R -e "devtools::install_local('monocle3.zip',upgrade = c('never'))"

RUN wget --tries=12 http://www.bioconductor.org/packages/release/bioc/src/contrib/GENIE3_1.18.0.tar.gz
RUN R CMD INSTALL GENIE3_1.18.0.tar.gz
RUN wget --tries=12 https://github.com/aertslab/SCopeLoomR/releases/download/v0.3.1/SCopeLoomR_0.3.1.tar.gz
RUN R CMD INSTALL SCopeLoomR_0.3.1.tar.gz
RUN R -e "install.packages('NMF')"
RUN wget --tries=12 https://github.com/aertslab/SCENIC/releases/download/v1.1.2/SCENIC_1.1.2.tar.gz
RUN R CMD INSTALL SCENIC_1.1.2.tar.gz

# 清理系统
RUN apt-get -y autoremove \
    && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# 入口 Google-Chrome -> http://127.0.0.1:8787