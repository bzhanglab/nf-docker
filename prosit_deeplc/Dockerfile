# Use an official Python runtime as a parent image
FROM python:3.10-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy the requirements file into the container
COPY prosit_deeplc_requirements.txt ./

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r prosit_deeplc_requirements.txt

# Command to run the Python script
CMD ["bash"]
