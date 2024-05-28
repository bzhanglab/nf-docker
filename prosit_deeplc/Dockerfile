# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Update package list and install procps
RUN apt-get update && apt-get install -y procps

# Clean up to reduce image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the requirements file into the container
COPY prosit_deeplc_requirements.txt ./

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r prosit_deeplc_requirements.txt

# Command to run the Python script
CMD ["bash"]
