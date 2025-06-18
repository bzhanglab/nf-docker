FROM openjdk:18-jdk-slim

# Update package list and install procps
RUN apt-get update && apt-get install -y procps

# Clean up to reduce image size
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the JAR file into the container
COPY pXg.v2.4.3.jar /app/pXg.jar
COPY inputconvertor.v1.1.0b.jar /app/inputconvertor.jar

CMD ["bash"]
