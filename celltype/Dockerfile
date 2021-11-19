FROM centos:7
RUN yum clean all && yum -y update
################## Java
RUN yum install -y java
##################  R (> version 3.4.3)
RUN yum install -y epel-release
RUN yum install -y R-core R-devel
#### utils
RUN yum install -y dos2unix wget mc
WORKDIR /opt
RUN cd /opt && wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh -O miniconda.sh \
    && chmod +x miniconda.sh && bash miniconda.sh -b -p miniconda
ENV PATH=/opt/miniconda/condabin:$PATH
RUN echo "conda activate" >> ~/.bashrc
RUN conda install -c conda-forge mamba
RUN conda install -y r-devtools
RUN source /opt/miniconda/bin/activate && Rscript -e "devtools::install_github('GfellerLab/EPIC')"
RUN mamba install -y -c bioconda -c conda-forge r-immunedeconv

