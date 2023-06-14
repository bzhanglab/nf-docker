FROM python:3.5
MAINTAINER Zhiao Shi <zhiao.shi@gmail.com>
LABEL software="pypgatk" \
    container="pypgatk" \
    software.version="0.0.19" \
    version="1"

RUN apt-get update
RUN apt-get install -y procps bedtools

WORKDIR /data
RUN mkdir -p /tool/source
COPY py-pgatk /tool/source/pypgatk
RUN cd /tool/source/pypgatk && pip install -r requirements.txt && python setup.py install

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV PATH=$PATH:/tool/source/pypgatk/
RUN chmod +x /tool/source/pypgatk/pypgatk/pypgatk_cli.py
