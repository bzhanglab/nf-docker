FROM ubuntu:18.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="Misc. scripts for rna-seq pipeline"
LABEL version="0.0.1"

RUN apt-get update && apt-get upgrade --yes && \
	apt-get install -y wget unzip curl build-essential
COPY src /usr/bin/

CMD ["/bin/bash"]
