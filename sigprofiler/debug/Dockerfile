FROM python:3.8

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

# Install the SigProfilerMatrixGenerator Python package in the virtual environment
RUN pip install --upgrade pip
RUN pip install SigProfilerPlotting==1.3.23
RUN pip install SigProfilerAssignment==0.1.6
RUN pip install SigProfilerMatrixGenerator==1.2.26
RUN pip install SigProfilerExtractor==1.1.24
RUN python -c "from SigProfilerMatrixGenerator import install as genInstall; genInstall.install('GRCh38', rsync=False, bash=True)"


RUN mkdir -p /opt/matplotlib_cache /opt/fontconfig_cache
# change perfmissions to allow R to write to these directories
RUN chmod 777 /opt/matplotlib_cache /opt/fontconfig_cache

# Set environment variables for matplotlib and fontconfig caches
ENV MPLCONFIGDIR=/opt/matplotlib_cache
ENV XDG_CACHE_HOME=/opt/fontconfig_cache

WORKDIR /usr/src/app
CMD ["/bin/bash"]