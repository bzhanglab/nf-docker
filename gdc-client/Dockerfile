FROM python:3.7.8
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

# Updating ubuntu and installing other necessary software
RUN apt-get clean all && rm -rf /var/lib/apt/lists/* \
	&& apt-get -y update --fix-missing \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
		tree \
		git \
		g++ \
		autoconf \
		automake \
		gcc \
		perl \
		libbz2-dev \
		liblzma-dev \
		libcurl4-gnutls-dev \
		libssl-dev \
		make \
		ncurses-dev \
		pkg-config \
		cmake  \
	&& apt-get clean

#install gdc-client
RUN wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip \
	&& unzip gdc-client_v1.6.1_Ubuntu_x64.zip \
	&& rm gdc-client_v1.6.1_Ubuntu_x64.zip \
	&& mv gdc-client /usr/bin/

#install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 \
	&& tar -xvjf htslib-1.15.1.tar.bz2 \
	&& rm htslib-1.15.1.tar.bz2 \
	&& cd htslib-1.15.1 \
	&& ./configure --prefix=/usr/local/ \
	&& make \
	&& make install \
	&& cd ..

#install Samtools
RUN wget https://github.com/samtools/samtools/archive/1.15.1.tar.gz \
	&& tar -xvzf 1.15.1.tar.gz \
	&& rm 1.15.1.tar.gz \
	&& cd samtools-1.15.1 \
	&& autoheader \
	&& autoconf -Wno-syntax \
	&& ./configure --prefix=/usr/local/ \
	&& make \
	&& make install \
	&& cd ..

#install VCFtools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz \
	&& tar -xvzf vcftools-0.1.16.tar.gz \
	&& rm vcftools-0.1.16.tar.gz \
	&& cd vcftools-0.1.16 \
	&& ./configure --prefix=/usr/local/ \
	&& make \
	&& make install \
	&& cd ..

#chmod of /opt/
RUN chmod -R 777 /opt/

#change working directory to /opt/
WORKDIR /opt/

#specify the command executed when the container is started
CMD ["/bin/bash"]