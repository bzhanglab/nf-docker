FROM openjdk:8
LABEL maintainer="Zhiao Shi<zhiao.shi@gmail.com>"
LABEL version="2.23.8"

# Install ant, git for building
RUN apt-get update && \
    apt-get --no-install-recommends install -y \
    git \
    r-base \
    pigz \
    ant && \
    apt-get clean autoclean && \
    apt-get autoremove -y 


RUN mkdir /usr/picard
COPY picard.jar /usr/picard/
