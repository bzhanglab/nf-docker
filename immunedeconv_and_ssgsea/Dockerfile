FROM rocker/tidyverse:4.3.0
LABEL maintainer="zhiao.shi@gmail.com"
WORKDIR /usr/src/app

# Install system dependencies for R packages
RUN apt-get update --fix-missing && apt-get install -y \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y aptitude && \
    aptitude install -y \
    libmagick++-dev \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*


RUN R -e "install.packages(c('corrplot', 'PerformanceAnalytics'), repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('GSVA')"
RUN R -e "BiocManager::install('gage')"
RUN R -e "BiocManager::install('progeny')"
RUN R -e "install.packages('estimate', repos='http://r-forge.r-project.org', dependencies=TRUE)"
RUN R -e "BiocManager::install('genefilter')"
RUN R -e "install.packages('purrr', repos='https://cran.r-project.org', dependencies=TRUE)"
RUN R -e "remotes::install_github('omnideconv/immunedeconv')"

CMD ["R"]
