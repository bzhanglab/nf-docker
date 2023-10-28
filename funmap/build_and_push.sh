#!/bin/bash

# Check if the version argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <version>"
    exit 1
fi

# Read the version from the command line argument
version="$1"

awk -v ver="$version" '{gsub(/funmap==.*/, "funmap=="ver)}1' Dockerfile > Dockerfile.tmp && mv Dockerfile.tmp Dockerfile

docker build -t registry.gitlab.com/bzhanglab/funmap:$version .
docker tag registry.gitlab.com/bzhanglab/funmap:$version registry.gitlab.com/bzhanglab/funmap:latest
docker push registry.gitlab.com/bzhanglab/funmap:$version
docker push registry.gitlab.com/bzhanglab/funmap:latest

