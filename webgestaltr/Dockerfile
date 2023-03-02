FROM rocker/tidyverse:4.0.2
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

ARG WGR_VER=0.4.5
RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -y curl
RUN R -e "install.packages(c('igraph','argparse', 'tools'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages(c('doRNG', 'doParallel', 'foreach', 'apcluster'), repos='https://cran.rstudio.com/')"
# svglite 2.1.1 does not compile
RUN R -e "install.packages(c('systemfonts'), repos='https://cran.rstudio.com/')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/svglite/svglite_2.0.0.tar.gz',repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/WebGestaltR_${WGR_VER}.tar.gz', repos=NULL, type='source')"
CMD ["R"]
