FROM jupyter/datascience-notebook:python-3.8.13
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

#specify the command executed when the container is started
RUN pip install lifelines==0.27.4
RUN pip install xgboost==1.6.2
RUN pip install statannotations==0.5.0
RUN pip install matplotlib==3.6.2
RUN pip install py4cytoscape==1.6.0

CMD ["/bin/bash"]
