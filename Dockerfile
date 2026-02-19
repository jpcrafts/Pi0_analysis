# syntax=docker/dockerfile:1

# Use a ROOT image based on Ubuntu 22.04 (jammy), not the oracular ones
FROM docker.io/rootproject/root:6.34.00-ubuntu24.04

ARG PI0_ANALYSIS_VERSION=dev
ARG VCS_REF=unknown
ARG BUILD_DATE=unknown

LABEL org.opencontainers.image.title="Pi0_analysis" \
    org.opencontainers.image.version="${PI0_ANALYSIS_VERSION}" \
    org.opencontainers.image.revision="${VCS_REF}" \
    org.opencontainers.image.created="${BUILD_DATE}"

# Update & install build tools + yaml-cpp
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
        build-essential \
        cmake \
        git \
        libyaml-cpp-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# Make sure CMake looks in /usr and ROOT’s own prefix
# (ROOT images usually already set CMAKE_PREFIX_PATH, but this is safe)
ENV CMAKE_PREFIX_PATH=/usr

WORKDIR /opt/Pi0_analysis

# Copy your project into the container
COPY . .

# Configure & build with RelWithDebInfo
RUN cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo && \
    cmake --build build --config RelWithDebInfo -j"$(nproc)"

# Put the build dir on PATH so all executables are easy to run
ENV PATH="/opt/Pi0_analysis/build:${PATH}"

# Default entrypoint: interactive shell
CMD ["/bin/bash"]
