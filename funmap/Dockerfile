FROM python:3.8
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

# Install the local package
RUN pip install funmap==0.1.18

ENV OMP_NUM_THREADS=1
ENV MKL_NUM_THREADS=1
#specify the command executed when the container is started
CMD ["/bin/bash"]
