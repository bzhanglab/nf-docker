FROM python:3.7

# Set the working directory inside the container
WORKDIR /app
# Copy the requirements file to the container
COPY pDeep3 .
RUN pip install -e .
CMD ["/bin/bash"]
