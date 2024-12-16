FROM python:2.7-slim
COPY . /usr/src/app
WORKDIR /usr/src/app

RUN mkdir -p /usr/share/man/man1
RUN apt-get update && apt-get install -y unzip vim wget bzip2 make gcc g++ zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev build-essential autoconf automake libtool openjdk-11-jre-headless procps

RUN wget https://github.com/arq5x/bedtools2/archive/master.zip
RUN unzip master.zip && rm master.zip

RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
RUN tar -xvf bwa-0.7.17.tar.bz2 && rm bwa-0.7.17.tar.bz2

RUN wget https://github.com/bioinfo-biols/CIRIquant/archive/refs/tags/v1.1.3.tar.gz
RUN tar -xvf v1.1.3.tar.gz && rm v1.1.3.tar.gz

RUN wget https://github.com/samtools/samtools/archive/1.9.zip
RUN unzip 1.9.zip && rm 1.9.zip

RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip
RUN unzip hisat2-2.1.0-Linux_x86_64.zip && rm hisat2-2.1.0-Linux_x86_64.zip

RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.3.tar.gz
RUN tar -xvf stringtie-2.0.3.tar.gz && rm stringtie-2.0.3.tar.gz
RUN wget https://github.com/samtools/htslib/archive/1.9.zip
RUN unzip 1.9.zip && rm 1.9.zip

# bedtools
WORKDIR /usr/src/app/bedtools2-master
RUN make && ln -s $PWD/bin/* /bin

# bwa
WORKDIR /usr/src/app/bwa-0.7.17
RUN make && ln -s $PWD/bwa /bin

# hisat
WORKDIR /usr/src/app/hisat2-2.1.0
RUN ln -s $PWD/hisat2-build $PWD/hisat2-build-{s,l} $PWD/hisat2-inspect-{s,l} $PWD/hisat2-align-{s,l} $PWD/hisat2 /bin

# stringtie
WORKDIR /usr/src/app/stringtie-2.0.3
RUN make clean release && ln -s $PWD/stringtie /bin

# htslib
WORKDIR /usr/src/app/htslib-1.9
RUN autoreconf && ./configure && make && make install

# samtools
WORKDIR /usr/src/app/samtools-1.9
RUN autoreconf && ./configure && make && make install
RUN ln -s $PWD/samtools /bin

# CIRIQuant
WORKDIR /usr/src/app/CIRIquant-1.1.3
RUN python setup.py install

# Create soft link link to Travis's location for things to work: 
RUN mkdir -p /home/travis/miniconda/envs/CIRIquant/ && ln -s /bin /home/travis/miniconda/envs/CIRIquant/


WORKDIR /usr/src/app
CMD ["/bin/bash"]
ENV PYTHON_EGG_CACHE=/tmp
