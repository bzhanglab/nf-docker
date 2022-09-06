#specifying the base image
FROM ubuntu:18.04
MAINTAINER Zhiao Shi zhiao.shi@gmail.com
# Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		git \
		g++ \
		make \
		ncurses-dev \
		pkg-config \
		wget \
		zlib1g-dev \
		cmake \
		bedtools \
	&& apt-get clean

#install STAR
RUN wget https://github.com/alexdobin/STAR/archive/2.7.10a.tar.gz \
	&& tar -xzf 2.7.10a.tar.gz \
	&& rm 2.7.10a.tar.gz \
	&& cd STAR-2.7.10a/source \
	&& make STAR \
	&& cd /

RUN wget https://github.com/williamritchie/IRFinder/archive/v1.3.1.tar.gz \
	&& tar -zxf v1.3.1.tar.gz \
	&& rm v1.3.1.tar.gz \
	&& chmod -R a+x IRFinder-1.3.1/bin

ENV PATH=$PATH:/STAR-2.7.10a/source
ENV PATH=$PATH:/IRFinder-1.3.1/bin

#change permission for the /IRFinder-1.2.3/ dir
RUN chmod -R 777 /IRFinder-1.3.1/

#change working directory to /opt/
WORKDIR /IRFinder-1.3.1/

#specify the command executed when the container is started
CMD ["/bin/bash"]
