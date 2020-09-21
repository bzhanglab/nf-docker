FROM python:3.7.8
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"

# Updating ubuntu and installing other necessary software
RUN apt-get clean all \
	&& apt-get -y update  \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		unzip \
		wget \
		build-essential \
		zlib1g-dev \
		libncurses5-dev \
		tree \
	&& apt-get clean

#install gdc-client
RUN wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.0_Ubuntu_x64-py3.7.zip \
	&& unzip gdc-client_v1.6.0_Ubuntu_x64-py3.7.zip \
	&& unzip gdc-client_v1.6.0_Ubuntu_x64.zip \
	&& rm gdc-client_v1.6.0_Ubuntu_x64.zip \
	&& rm gdc-client_v1.6.0_Ubuntu_x64-py3.7.zip \
	&& mv gdc-client /usr/bin/

#specify the command executed when the container is started
CMD ["/bin/bash"]