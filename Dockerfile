# Minimal Docker image for ViralMSA using Ubuntu base
FROM ubuntu:20.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake g++ g++-10 gcc-10 libboost-all-dev libgsl-dev libjemalloc-dev make unzip wget zlib1g-dev && \

    # Install Bowtie2 v2.4.3
    wget "https://github.com/BenLangmead/bowtie2/releases/download/v2.4.3/bowtie2-2.4.3-source.zip" && \
    unzip bowtie2-*-source.zip && \
    cd bowtie2-* && \
    make && \
    make install && \
    cd .. && \
    rm -rf bowtie2-* && \

    # Install DRAGMAP v1.2.1
    wget -qO- "https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.2.1.tar.gz" | tar -zx && \
    cd DRAGMAP-* && \
    HAS_GTEST=0 make CFLAGS:= && \
    HAS_GTEST=0 make install && \
    cd .. && \
    rm -rf DRAGMAP-* && \

    # Install HISAT2 v2.2.1
    wget -qO- "https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz" | tar -zx && \
    cd hisat2-* && \
    make && \
    mv hisat2 hisat2-* hisat2_*.py /usr/local/bin/ && \
    cd .. && \
    rm -rf hisat2-* && \

    # Install Minimap2 v2.24
    wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.24.tar.gz" | tar -zx && \
    cd minimap2-* && \
    make && \
    chmod a+x minimap2 && \
    mv minimap2 /usr/local/bin/minimap2 && \
    cd .. && \
    rm -rf minimap2-* && \

    # Install NGMLR v0.2.7
    wget -qO- "https://github.com/philres/ngmlr/archive/refs/tags/v0.2.7.tar.gz" | tar -zx && \
    cd ngmlr-* && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make && \
    mv ../bin/ngmlr-*/ngmlr /usr/local/bin/ngmlr && \
    cd ../.. && \
    rm -rf ngmlr-* && \

    # Install STAR v2.7.5c
    wget -qO- "https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz" | tar -zx && \
    mv STAR-*/bin/Linux_*_static/* /usr/local/bin/ && \
    rm -rf STAR-* && \

    # Install Unimap (latest)
    wget "https://github.com/lh3/unimap/archive/refs/heads/master.zip" && \
    unzip master.zip && \
    cd unimap-master && \
    make && \
    mv unimap /usr/local/bin/unimap && \
    cd .. && \
    rm -rf master.zip unimap-master && \

    # Install wfmash v0.7.0
    wget -qO- "https://github.com/ekg/wfmash/releases/download/v0.7.0/wfmash-v0.7.0.tar.gz" | tar -zx && \
    cd wfmash-* && \
    cmake -H. -Bbuild -DCMAKE_C_COMPILER="$(which gcc-10)" -DCMAKE_CXX_COMPILER="$(which g++-10)" && \
    cmake --build build -- && \
    mv build/bin/wfmash /usr/local/bin/wfmash && \
    cd .. && \
    rm -rf wfmash-* && \

    # Install Windowmap v2.03
    wget -qO- "https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz"  | tar -zx && \
    cd Winnowmap-* && \
    make && \
    mv bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf Winnowmap-* && \

    # Set up ViralMSA
    wget -O /usr/local/bin/ViralMSA.py "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py" && chmod a+x /usr/local/bin/ViralMSA.py && \
    
    # Clean up
    rm -rf /root/.cache /tmp/*
