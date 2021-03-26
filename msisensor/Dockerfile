FROM ubuntu:18.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="msisensor (https://github.com/ding-lab/msisensor)"
LABEL version="0.6"

RUN apt-get update \
    && apt-get upgrade --yes \
	&& apt-get install -y --no-install-recommends wget build-essential curl ca-certificates zlib1g-dev  ncurses-dev\
	&& rm -rf /var/lib/apt/lists/*

RUN mkdir /build && cd /build && curl -OL https://github.com/ding-lab/msisensor/archive/0.6.tar.gz && tar xzf 0.6.tar.gz \
    && cd msisensor-0.6 && make && cp msisensor /usr/local/bin/msisensor

CMD ["/bin/bash"]
