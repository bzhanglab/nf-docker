FROM rocker/tidyverse:4.0.2
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL version="1.0.0"

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -y curl

RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install(version = "3.12", ask = FALSE); BiocManager::install("maftools")' \
    && Rscript -e "stopifnot(all(c('BiocManager', 'maftools') %in% installed.packages()[, 'Package']))"
COPY src/* /usr/bin/
CMD ["/bin/bash"]
