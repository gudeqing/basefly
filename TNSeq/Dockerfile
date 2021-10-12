FROM ensemblorg/ensembl-vep:release_104.3
ARG VERSION=202010.02
USER root
# Download the software from the permalink
RUN apt-get update && apt-get install -y wget python3.8 nano python3.8-dev libgraphviz-dev --fix-missing && \
    mkdir -p /opt/sentieon/ && \
    wget -nv -O - "https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-${VERSION}.tar.gz" | \
      tar -zxf - -C /opt/sentieon/

CMD ["/bin/bash"]

ENV PATH /opt/sentieon/sentieon-genomics-${VERSION}/bin/:$PATH
ENV SENTIEON_INSTALL_DIR=/opt/sentieon-genomics-${VERSION} \
    SENTIEON_LICENSE=sentieon-lic.enigma2:8990

# python packages
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN apt install -y python3-pip
RUN pip3 install -i $PIP_SOURCE dataclasses typing_extensions scipy numpy pandas biopython psutil pygraphviz
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl
