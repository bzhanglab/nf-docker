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
RUN R -e "BiocManager::install('org.Hs.eg.db', version = '${BIOC_VER}')"
RUN R -e "install.packages(c('doParallel','foreach', 'R2HTML'), repos='https://cran.rstudio.com/')"
RUN R -e "BiocManager::install('biomaRt', version = '${BIOC_VER}')"
RUN R -e "install.packages(c('seriation'), repos='https://cran.rstudio.com/')"
RUN R -e "BiocManager::install('impute', version = '${BIOC_VER}')"
RUN R -e "BiocManager::install('preprocessCore', version = '${BIOC_VER}')"
RUN R -e "install.packages(c('WGCNA'), repos='https://cran.rstudio.com/')"
COPY NetSAM_1.38.1.tar.gz /tmp/
RUN R CMD INSTALL /tmp/NetSAM_1.38.1.tar.gz
CMD ["R"]
