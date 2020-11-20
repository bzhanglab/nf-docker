FROM jupyter/datascience-notebook:python-3.7.6
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL version=1.0.1

#specify the command executed when the container is started
USER root
RUN apt-get update -y
RUN apt-get install -y tcsh parallel
RUN pip install biopython
COPY binding_prediction.py /usr/local/bin
CMD ["/bin/bash"]
