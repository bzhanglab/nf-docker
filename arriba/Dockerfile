FROM ubuntu:bionic
LABEL maintainer="zhiao.shi@gmail.com"

# install dependencies
RUN export DEBIAN_FRONTEND=noninteractive && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends build-essential samtools r-base wget ca-certificates libcurl4-openssl-dev libxml2-dev vim

# install version of STAR that supports --chimMultimapNmax and --chimOutType WithinBAM
RUN wget -q -O - 'https://github.com/alexdobin/STAR/archive/2.7.6a.tar.gz' | tar -xzf - && \
    cp -p STAR-2.7.6a/bin/Linux_x86_64/STAR /usr/local/bin/

# install arriba
RUN wget -q -O - "https://github.com/suhrig/arriba/releases/download/v2.0.0/arriba_v2.0.0.tar.gz" | tar -xzf -

# this scripts is updated base on the one comes with arriba
# add GRCh38.p13 and GENCODEV34_basic
COPY download_references.sh /usr/local/bin 
RUN chmod a+x /usr/local/bin/download_references.sh

# modified run script 
COPY run_arriba.sh /arriba_v2.0.0
RUN chmod a+x /arriba_v2.0.0/run_arriba.sh
RUN ln -s /arriba_v2.0.0/run_arriba.sh /usr/local/bin
RUN ln -s /arriba_v2.0.0 /arriba
RUN ln -s /arriba_v2.0.0/arriba /usr/local/bin