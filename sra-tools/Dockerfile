FROM ubuntu:20.04 AS builder
RUN apt-get update && apt-get install -y build-essential zlib1g-dev curl && \
    curl -L https://github.com/madler/pigz/archive/refs/tags/v2.8.tar.gz -o pigz.tar.gz && \
    tar -xzf pigz.tar.gz && \
    cd pigz-2.8 && \
    make

# Stage 2: Final image
FROM quay.io/biocontainers/sra-tools:3.2.1--h4304569_1

USER root
# Copy zlib shared library from builder
COPY --from=builder /lib/x86_64-linux-gnu/libz.so.1* /usr/lib/

# Copy compiled pigz binary from builder
COPY --from=builder /pigz-2.8/pigz /usr/local/bin/pigz
RUN chmod +x /usr/local/bin/pigz
