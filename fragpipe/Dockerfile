FROM ubuntu:20.04
LABEL maintainer="zhiao.shi@gmail.com"
LABEL version="18.0"
LABEL about="Computational tools enabling comprehensive analysis of mass spectrometry-based proteomics data"

USER root
RUN apt-get update && apt-get -y install wget unzip default-jre-headless && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN ZIP=FragPipe-18.0.zip && wget https://github.com/Nesvilab/FragPipe/releases/download/18.0/${ZIP} -O /tmp/FragPipe-18.0.zip &&  unzip /tmp/${ZIP} -d /usr/local

ENV PATH=/usr/local/fragpipe/bin:${PATH}
