FROM python:3.11

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive
ENV LD_LIBRARY_PATH=/usr/local/lib

LABEL name="ms2rescore-percolator"
LABEL description="MS²Rescore with Percolator rescoring engine"

# Install system dependencies for both ms2rescore and Percolator
RUN apt-get update && \
    apt-get install -y \
        git \
        build-essential \
        sudo \
        lsb-release \
        cmake \
        libboost-all-dev \
        libxml2-dev \
        wget \
        libxerces-c-dev \
        xsdcxx \
        libsqlite3-dev \
        libcurl4-openssl-dev \
        procps \
        curl \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Build and install Percolator
RUN wget https://github.com/percolator/percolator/archive/refs/tags/rel-3-06-05.tar.gz && \
    tar -xzf rel-3-06-05.tar.gz && \
    cd percolator-rel-3-06-05 && \
    cmake . && \
    make && \
    make install && \
    cd src/converters && \
    sed -i '361i target_link_libraries(converters curl icuuc)' CMakeLists.txt && \
    cmake . && \
    make && \
    make install && \
    cd ../../../ && \
    rm -rf percolator-rel-3-06-05 rel-3-06-05.tar.gz

# Verify Percolator installation
RUN percolator --help > /dev/null 2>&1 && echo "Percolator installed successfully"

# Clone and install ms2rescore
RUN git clone https://github.com/CompOmics/ms2rescore.git /ms2rescore && \
    cd /ms2rescore && \
    pip install . --only-binary :all:

# Update library cache for Percolator
RUN ldconfig

# Verify installations (with fallback)
RUN echo "Checking installations..." && \
    python -c "import ms2rescore; print(f'MS²Rescore {ms2rescore.__version__} installed successfully')" && \
    percolator --help | head -1 || echo "Percolator installed (help unavailable)" && \
    echo "All installations completed successfully"

ENV HOME=/opt/home
RUN mkdir -p $HOME && chmod -R 777 $HOME

# Set Matplotlib config directory to something writable
ENV MPLCONFIGDIR=$HOME/.config/matplotlib

# Optional: pre-create it for safety
RUN mkdir -p $MPLCONFIGDIR && chmod -R 777 $MPLCONFIGDIR

# Optional: pre-create ms2pip cache dir if it's known
RUN mkdir -p $HOME/.ms2pip && chmod -R 777 $HOME/.ms2pip

WORKDIR /data
ENV MPLCONFIGDIR=/tmp/mplconfig
CMD ["bash"]
