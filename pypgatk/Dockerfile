FROM python:3.9
# FROM python:3.7
LABEL MAINTAINER="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL software="pypgatk" \
    container="pypgatk" \
    software.version="0.0.23" \
    version="1"

RUN apt-get update
RUN apt-get install -y procps bedtools 
RUN pip install jupyter

# WORKDIR /data
# RUN mkdir -p /tool/source
# COPY py-pgatk /tool/source/pypgatk
# RUN cd /tool/source/pypgatk && pip install -r requirements.txt && python setup.py install

# ENV LC_ALL=C.UTF-8
# ENV LANG=C.UTF-8
# ENV PATH=$PATH:/tool/source/pypgatk/
# RUN chmod +x /tool/source/pypgatk/pypgatk/pypgatk_cli.py

WORKDIR /data
RUN mkdir -p /tool/source && cd /tool/source

RUN git config --global http.sslVerify false
RUN git clone --depth 1 https://github.com/bigbio/py-pgatk.git /tool/source
RUN cd /tool/source/ && pip3 install -r requirements.txt && python3 setup.py install

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH=$PATH:/tool/source/pypgatk/
RUN chmod +x /tool/source/pypgatk/pypgatk_cli.py