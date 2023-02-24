FROM alpine:3.17

RUN apk add --no-cache gcc musl-dev glib-dev make
# copy the source code to the container, do i need to do this recursively?
COPY ice-1.0.2.tgz /tmp/
# build and install the source code
RUN cd /tmp && tar xvfz ice-1.0.2.tgz && cd ice && ./configure && make && make install