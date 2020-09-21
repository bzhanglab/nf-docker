FROM ubuntu:18.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="RSEM: accurate quantification of gene and isoform expression from RNA-Seq data"
LABEL version="1.3.3"

#updating ubuntu and installing other necessary software
RUN apt-get update --fix-missing \
	&& apt-get install -y \
		wget \
		unzip \
		bzip2 \
		ca-certificates \
		build-essential \
		checkinstall \
		zlib1g-dev \
		libncurses5-dev \
		python3

RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.1/bowtie2-2.4.1-linux-x86_64.zip \
	&& unzip bowtie2-2.4.1-linux-x86_64.zip \
	&& rm bowtie2-2.4.1-linux-x86_64.zip \
	&& mv bowtie2-2.4.1-linux-x86_64 /opt/ \ 
	&& ln -s /opt/bowtie2-2.4.1-linux-x86_64 /opt/bowtie2

#install RSEM
RUN wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz \
	&& tar xzvf v1.3.3.tar.gz \
	&& rm v1.3.3.tar.gz \
	&& cd RSEM-1.3.3 \
	&& make \
	&& cd .. \
	&& mv RSEM-1.3.3 /opt/ \
  && ln -s /opt/RSEM-1.3.3 /opt/RSEM

# http://bugs.python.org/issue19846
ENV LANG C.UTF-8

CMD ["/bin/bash"]
