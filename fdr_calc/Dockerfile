FROM proteomics/pga:latest
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -y curl
COPY fdr_calc.R /usr/src 
CMD ["/bin/bash"]
