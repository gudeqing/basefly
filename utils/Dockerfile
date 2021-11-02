FROM python:3.9-slim-bullseye
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN pip install -i $PIP_SOURCE multiqc plotly colorlover fastcluster pandas_profiling openpyxl
RUN pip install -i $PIP_SOURCE bokeh
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl


