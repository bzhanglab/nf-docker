FROM r-base:4.4.0

# Install necessary system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and required packages
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"
RUN R -e "BiocManager::install(version = '3.19')"
RUN R -e "BiocManager::install('maftools')"

# Install CRAN packages
RUN R -e "install.packages(c('tidyr', 'dplyr'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('R.utils', repos='http://cran.rstudio.com/')"

# Set the working directory
WORKDIR /usr/src/app

# Copy your R script or project files into the container
#COPY . .

# Default command to run when starting the container
CMD ["R"]
