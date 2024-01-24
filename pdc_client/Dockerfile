# Use Ubuntu 20.04 as the base image
FROM ubuntu:20.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install Python 3.8
RUN apt-get update && apt-get install -y python3.8 python3-pip

# Set the working directory in the container
WORKDIR /app

# Copy the pdc_client executable to the container
COPY pdc-client /usr/bin/pdc-client

# Give execute permission to the pdc_client
RUN chmod +x /usr/bin/pdc-client

# Command to run when the container starts
CMD ["/bin/bash"]

