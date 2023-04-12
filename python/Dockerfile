FROM jupyter/datascience-notebook:python-3.7.6
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

#specify the command executed when the container is started
RUN pip install lifelines==0.27.4
RUN pip install xgboost
RUN pip install statannotations

CMD ["/bin/bash"]
