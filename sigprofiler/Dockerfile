# Use the R base image with R 4.4 installed
FROM r-base:4.4.0

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    wget \
    zlib1g-dev \
    libncurses5-dev \
    libgdbm-dev \
    libnss3-dev \
    libsqlite3-dev \
    libreadline-dev \
    libffi-dev \
    libbz2-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Python from source with --enable-shared
RUN wget https://www.python.org/ftp/python/3.10.4/Python-3.10.4.tgz && \
    tar -xf Python-3.10.4.tgz && \
    cd Python-3.10.4 && \
    ./configure --enable-shared && \
    make && \
    make install && \
    ldconfig && \
    cd .. && \
    rm -rf Python-3.10.4 Python-3.10.4.tgz

# Create and activate a virtual environment for Python
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Install the SigProfilerMatrixGenerator Python package in the virtual environment
RUN pip install --upgrade pip
RUN pip install SigProfilerMatrixGenerator

# Install R packages
RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reticulate', repos='http://cran.rstudio.com/')"

# Configure reticulate to use the created Python virtual environment
RUN R -e "library(reticulate); use_virtualenv('/opt/venv', required = TRUE)"

# Install SigProfilerMatrixGeneratorR from GitHub
RUN R -e "devtools::install_github('AlexandrovLab/SigProfilerMatrixGeneratorR')"

# Install the reference genome for SigProfilerMatrixGeneratorR
RUN R -e "library('SigProfilerMatrixGeneratorR'); install('GRCh38', rsync=FALSE, bash=TRUE)"
RUN R -e "devtools::install_github('AlexandrovLab/SigProfilerPlottingR')"
RUN R -e "devtools::install_github('AlexandrovLab/SigProfilerExtractorR')"

ENV RETICULATE_PYTHON="/opt/venv/bin/python"

WORKDIR /usr/src/app
CMD ["R"]
