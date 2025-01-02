FROM ubuntu:18.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="CircRNA Identifier: A de novo circular RNA identification tool"
LABEL version="2.0.6"

RUN apt-get update && apt-get upgrade --yes && \
	apt-get install -y wget unzip curl build-essential zlib1g-dev
RUN mkdir /usr/src/ciri && \
	cd /usr/src/ciri && \
	wget https://sourceforge.net/projects/ciri/files/CIRI2/CIRI_v2.0.6.zip && \
	unzip CIRI_v2.0.6.zip && \
	cp CIRI_v2.0.6/CIRI2.pl /usr/src && \
	cd /usr/src && \
	wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
	tar -xjf bwa-0.7.17.tar.bz2 && \
	cd bwa-0.7.17 && \
	make && \
	cp bwa /usr/local/bin/

CMD ["/bin/bash"]
