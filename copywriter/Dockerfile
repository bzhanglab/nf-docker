# dockerfile for bzhanglab:copywriter 
FROM rocker/tidyverse:latest
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="R package CopywriteR for detection of copy number variants and aberrations"
RUN apt-get update -qq && apt-get install -y --no-install-recommends \
    libbz2-dev \
    liblzma-dev
RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(version = "3.11", ask = FALSE); BiocManager::install("CopywriteR", version = "3.11")' \
    && Rscript -e "stopifnot(all(c('BiocManager', 'CopywriteR') %in% installed.packages()[, 'Package']))"
CMD ["R"]
