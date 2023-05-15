FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        ca-certificates \
        language-pack-en \
        locales \
        locales-all \
        wget \
        unzip 

RUN wget https://github.com/lmsac/phosphoRS-cli/releases/download/v1.0.0/phosphoRS-cli.zip \
    # unzip to a folder named phosphoRS-cli
    && mkdir phosphoRS-cli \
    && mv phosphoRS-cli.zip phosphoRS-cli \
    && cd phosphoRS-cli \
    && unzip phosphoRS-cli.zip \
    && rm phosphoRS-cli.zip \
    && cd .. \
    && mv phosphoRS-cli /usr/local/bin/
# Install Wine
RUN dpkg --add-architecture i386 && \
    mkdir -pm755 /etc/apt/keyrings && \
    wget -O /etc/apt/keyrings/winehq-archive.key https://dl.winehq.org/wine-builds/winehq.key && \
    wget -nc -P /etc/apt/sources.list.d/ https://dl.winehq.org/wine-builds/ubuntu/dists/$(lsb_release -sc)/winehq-$(lsb_release -sc).sources && \
    apt-get update -y && \
    apt-get install -y --install-recommends winehq-stable=8.0.1~focal-1