FROM rocker/tidyverse:3.6.3
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

RUN apt-get update && apt-get upgrade --yes && \
    apt-get install -y curl && apt-get install -y libmagick++-dev
RUN R -e "devtools::install_version('reshape2', version='1.4.4', upgrade='never', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('gridExtra', version='2.3', upgrade='never', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "devtools::install_version('kableExtra', version='1.3.1', upgrade='never', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')"
RUN R -e "devtools::install_version('mutSignatures', version='2.1.1', upgrade='never', dependencies=TRUE, repos='http://cran.rstudio.com/')"

CMD ["R"]