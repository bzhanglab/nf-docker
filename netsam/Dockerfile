FROM rocker/tidyverse:4.2
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

ARG BIOC_VER=3.16
RUN apt-get update && apt-get install -y \
    curl \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev
RUN R -e "install.packages(c('igraph','argparse', 'tools', 'DBI'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('GO.db', version = '${BIOC_VER}')"
# netsam version is 1.38.0
RUN R -e "BiocManager::install('NetSAM', version = '${BIOC_VER}')"
CMD ["R"]
