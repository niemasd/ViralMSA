# Minimal Docker image for ViralMSA using Ubuntu base
FROM ubuntu:20.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y g++ make unzip wget zlib1g-dev

# Install Bowtie2 v2.4.3
RUN wget "https://github.com/BenLangmead/bowtie2/releases/download/v2.4.3/bowtie2-2.4.3-source.zip" && \
    unzip bowtie2-2.4.3-source.zip && \
    cd bowtie2-2.4.3 && \
    make && \
    make install && \
    cd .. && \
    rm -rf bowtie2-2.4.3 bowtie2-2.4.3-source.zip

# Install HISAT2 v2.2.1
RUN wget -qO- "https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz" | tar -zx && \
    cd hisat2-* && \
    make && \
    mv hisat2 hisat2-* hisat2_*.py /usr/local/bin/ && \
    cd .. && \
    rm -rf hisat2-*

# Install Minimap2 v2.24
RUN wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz" | tar -zx && \
    cd minimap2-* && \
    make && \
    chmod a+x minimap2 && \
    mv minimap2 /usr/local/bin/minimap2 && \
    cd .. && \
    rm -rf minimap2-*

# Install STAR v2.7.5c
RUN wget -qO- "https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz" | tar -zx && \
    mv STAR-*/bin/Linux_*_static/* /usr/local/bin/ && \
    rm -rf STAR-*

# Install Unimap (latest)
RUN wget "https://github.com/lh3/unimap/archive/refs/heads/master.zip" && \
    unzip master.zip && \
    cd unimap-master && \
    make && \
    mv unimap /usr/local/bin/unimap && \
    cd .. && \
    rm -rf master.zip unimap-master

# Set up ViralMSA
RUN wget -O /usr/local/bin/ViralMSA.py "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py" && chmod a+x /usr/local/bin/ViralMSA.py

# Clean up
RUN rm -rf /root/.cache && \
    rm -rf /tmp/*
