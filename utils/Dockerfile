FROM python:3.9-slim-bullseye
ENV PIP_SOURCE https://mirrors.aliyun.com/pypi/simple
RUN pip install -i $PIP_SOURCE multiqc plotly colorlover fastcluster pandas_profiling openpyxl
RUN pip install -i $PIP_SOURCE bokeh pysam
RUN pip3 install https://github.com/gudeqing/xcmds/raw/master/dist/xcmds-1.5.0-py3-none-any.whl
ADD clusterheatmap.py /bin/clusterheatmap.py
RUN pip install ipython
RUN mkdir -p ~/.ipython/profile_default/startup/ && \
    echo " import os; import scipy.stats as stats; \
    import matplotlib; matplotlib.use('agg'); import pandas as pd; import numpy as np; \
    from matplotlib import pyplot as plt; import seaborn as sns" > \
    ~/.ipython/profile_default/startup/00-init_imports.py
RUN apt-get update && apt-get install -y vim
