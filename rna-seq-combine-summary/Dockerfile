FROM rocker/tidyverse:latest
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL version="1.0.0"

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -y curl
RUN install2.r --error corrplot purrr factoextra BiocManager
RUN Rscript -e 'BiocManager::install(c("edgeR", "preprocessCore"))'
COPY src/combine_and_summary.r /usr/src

CMD ["R"]
