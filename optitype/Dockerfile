FROM biocontainers/biocontainers:latest
LABEL maintaainer="Zhiao Shi<zhiao.shi@gmail.com>"

USER root
# install
RUN apt-get update && apt-get install -y software-properties-common \
&& apt-get update && apt-get install -y \
    gcc-4.9 \
    g++-4.9 \
    coinor-cbc \
    zlib1g-dev \
    libbz2-dev \
&& update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.9 \
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean \
&& apt-get purge

#HLA Typing
#OptiType dependecies
RUN curl -O https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.8/hdf5-1.8.21/bin/hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz \
    && tar -xvf hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz \
    && mv hdf5/bin/* /usr/local/bin/ \
    && mv hdf5/lib/* /usr/local/lib/ \
    && mv hdf5/include/* /usr/local/include/ \
    && mv hdf5/share/* /usr/local/share/ \
    && rm -rf hdf5 \
    && rm -f hdf5-1.8.21-Std-centos7-x86_64-shared_64.tar.gz

ENV LD_LIBRARY_PATH /usr/local/lib:$LD_LIBRARY_PATH
ENV HDF5_DIR /usr/local/

RUN pip install --upgrade pip && pip install \
    numpy \
    pyomo \
    pysam \
    matplotlib \
    tables \
    pandas \
    future

#installing razers3
RUN git clone https://github.com/seqan/seqan.git seqan-src \
        && cd seqan-src \
        && cmake -DCMAKE_BUILD_TYPE=Release \
        && make razers3 \
        && mv bin/razers3 /opt/ \
        && cd .. \
        && rm -rf seqan-src

#installing optitype form git repository (version Dec 09 2015) and wirtig config.ini
RUN git clone https://github.com/FRED-2/OptiType.git \
        && mv OptiType/ /opt/ \
        && chmod 777 /opt/OptiType/OptiTypePipeline.py \
        && echo "[mapping]\n\
razers3=/opt/razers3 \n\
threads=1 \n\
\n\
[ilp]\n\
solver=cbc \n\
threads=1 \n\
\n\
[behavior]\n\
deletebam=true \n\
unpaired_weight=0 \n\
use_discordant=false\n" >> /opt/OptiType/config.ini
ENV PATH=/opt/OptiType:$PATH
#chmod of /home/user and change working directory
RUN mkdir /.config
RUN chmod -R 777 /opt/
RUN chmod -R 700 /.config
RUN chmod 777 /home/biodocker
#change working directory
WORKDIR /opt/
