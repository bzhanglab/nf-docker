FROM ubuntu:18.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="Burrow-Wheeler Aligner for short-read alignment"
LABEL version="0.7.17"

RUN apt-get update \
    && apt-get upgrade --yes \
	&& apt-get install -y --no-install-recommends wget build-essential zlib1g-dev \
	&& rm -rf /var/lib/apt/lists/*
RUN mkdir /usr/bwa \
	&& cd /usr/bwa \
    && wget --no-check-certificate https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
	&& tar xvjf bwa-0.7.17.tar.bz2 \
	&& rm -f bwa-0.7.17.tar.bz2 \
	&& cd bwa-0.7.17 \
	&& make \
	&& cp bwa /usr/local/bin  

CMD ["/bin/bash"]