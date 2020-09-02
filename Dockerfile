# ViralMSA minimal Docker image using
FROM python:3.8-slim-buster
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y bzip2 curl && \
    pip3 install biopython

# Install Minimap2 (2.17)
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2 | tar xj && \
    mv minimap2-*/minimap2 /usr/local/bin && rm -rf minimap2-*
