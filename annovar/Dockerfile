FROM ubuntu:20.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

#specify the command executed when the container is started

# copy annovar.tgz to the container
COPY annovar.tgz /usr/src/annovar.tgz

ENV DEBIAN_FRONTEND=noninteractive

# Update the package list and install required packages
RUN apt-get update && apt-get install -y \
    perl \
    cpanminus \
    make \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Upgrade Pod::Simple to the required version
RUN cpanm Pod::Simple@3.40

# Install the Pod::Usage Perl module
RUN cpanm Pod::Usage

# set the working directory
#WORKDIR /usr/src/app/annovar

CMD ["/bin/bash"]
