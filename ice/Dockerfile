FROM ubuntu:20.04

ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y gcc musl-dev libglib2.0-dev make wget

# copy the source code to the container, do i need to do this recursively?
COPY ice-1.0.2.tar.gz /tmp/
# build and install the source code
RUN cd /tmp && tar xvfz ice-1.0.2.tar.gz && cd ice && ./configure && make && make install
