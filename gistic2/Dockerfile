FROM ubuntu:14.04
LABEL maintainer="Zhiao Shi <zhiao.shi@gmail.com>"
LABEL description="Genomic Identification of Significant Targets in Cancer"
LABEL version="2.0.23"

# Note: FROM java and FROM r-base work too but take much longer apt-get updating

RUN apt-get update && apt-get upgrade --yes && \ 
    apt-get install -y wget && \
    apt-get install --yes bc vim libxpm4 libXext6 libXt6 libXmu6 libXp6 zip unzip  build-essential

RUN mkdir /home/gistic
WORKDIR /home/gistic
RUN wget ftp://ftp.broadinstitute.org/pub/GISTIC2.0/GISTIC_2_0_23.tar.gz \
    && tar xvf GISTIC_2_0_23.tar.gz
RUN mkdir /home/gistic/MCRInstaller
RUN cd /home/gistic/MCRInstaller \
    && wget https://www.mathworks.com/supportfiles/downloads/R2014a/deployment_files/R2014a/installers/glnxa64/MCR_R2014a_glnxa64_installer.zip \
    && unzip MCR_R2014a_glnxa64_installer.zip
RUN cd MCRInstaller \
    && /home/gistic/MCRInstaller/install -mode silent -agreeToLicense yes  -destinationFolder /home/gistic/MATLAB_Compiler_Runtime
RUN chmod o+rx /home/gistic
ENV LD_LIBRARY_PATH /home/gistic/MATLAB_Compiler_Runtime/v83/runtime/glnxa64:/home/gistic/MATLAB_Compiler_Runtime/v83/bin/glnxa64:/home/gistic/MATLAB_Compiler_Runtime/v83/sys/os/glnxa64:${LD_LIBRARY_PATH}
ENV XAPPLRESDIR /home/gistic/MATLAB_Compiler_Runtime/v83/X11/app-defaults
ENV PATH=/home/gistic:${PATH}
COPY gencode.v34.basic.annotation.mat /home/gistic/refgenefiles
CMD ["/bin/bash"]
