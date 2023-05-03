FROM python:3.8
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

ENV VERSION=0.1.10
COPY funmap /usr/src/app/funmap

# Set the working directory to the package directory
WORKDIR /usr/src/app/funmap

# Install the local package
RUN pip install .

ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
#specify the command executed when the container is started
CMD ["/bin/bash"]
