FROM rocker/tidyverse:4.5.3@sha256:3a4f24ca1159bd822383cfb85f39d5602d9e2f447ace9271ab1330a207249328
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL version="1.0.3"

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -y curl
RUN install2.r --error corrplot purrr factoextra BiocManager
RUN Rscript -e 'BiocManager::install(c("edgeR", "preprocessCore"))'
COPY src/combine_and_summary.r /usr/src
COPY src/combine_and_summary_polya.r /usr/src

CMD ["R"]
