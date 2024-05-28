FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y git build-essential sudo lsb-release cmake \
    libboost-all-dev libxml2-dev wget libxerces-c-dev xsdcxx libsqlite3-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/percolator/percolator/archive/refs/tags/rel-3-06-05.tar.gz && \
    tar -xzf rel-3-06-05.tar.gz && \
    cd percolator-rel-3-06-05 && \
    cmake . && \
    make && \
    make install  
RUN cd percolator-rel-3-06-05/src/converters && \
    sed -i '361i target_link_libraries(converters curl icuuc)' CMakeLists.txt && \
    cmake . && \
    make && \
    make install

CMD ["/bin/bash"]
