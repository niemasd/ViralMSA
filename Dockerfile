# Minimal Docker image for ViralMSA using Ubuntu base
FROM ubuntu:20.04
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get install -y cmake g++ g++-10 gcc-10 git libboost-all-dev libbz2-dev libcurl4-openssl-dev libgsl-dev libjemalloc-dev liblzma-dev make pkg-config python3 python3-pip unzip wget zlib1g-dev && \

    # Install Python packages
    pip install --no-cache-dir biopython && \

    # Install htslib
    wget -qO- "https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2" | tar -xj && \
    cd htslib-* && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-* && \

    # Install Bowtie2
    wget -qO- "https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.3.tar.gz" | tar -zx && \
    cd bowtie2-* && \
    make && \
    make install && \
    cd .. && \
    rm -rf bowtie2-* && \

    # Install BWA
    wget -qO- "https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2" | tar -xj && \
    cd bwa-* && \
    sed -i 's/const uint8_t rle_auxtab\[8\];/\/\/const uint8_t rle_auxtab\[8\];/g' rle.h && \
    make && \
    mv bwa /usr/local/bin/bwa && \
    cd .. && \
    rm -rf bwa-* && \

    # Install DRAGMAP
    wget -qO- "https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.3.0.tar.gz" | tar -zx && \
    cd DRAGMAP-* && \
    HAS_GTEST=0 make CFLAGS:= && \
    HAS_GTEST=0 make install && \
    cd .. && \
    rm -rf DRAGMAP-* && \

    # Install HISAT2
    wget -qO- "https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.1.tar.gz" | tar -zx && \
    cd hisat2-* && \
    make && \
    mv hisat2 hisat2-* hisat2_*.py /usr/local/bin/ && \
    cd .. && \
    rm -rf hisat2-* && \

    # Install Minimap2
    wget -qO- "https://github.com/lh3/minimap2/archive/refs/tags/v2.28.tar.gz" | tar -zx && \
    cd minimap2-* && \
    make && \
    chmod a+x minimap2 && \
    mv minimap2 /usr/local/bin/minimap2 && \
    cd .. && \
    rm -rf minimap2-* && \

    # Install mm2-fast
    wget -qO- "https://github.com/bwa-mem2/mm2-fast/releases/download/mm2-fast-v1.0/Source_code_including_submodules.tar.gz" | tar -zx && \
    cd mm2-fast* && \
    make && \
    mv minimap2 /usr/local/bin/mm2-fast && \
    cd .. && \
    rm -rf mm2-fast* && \

    # Install NGMLR
    wget -qO- "https://github.com/philres/ngmlr/archive/refs/tags/v0.2.7.tar.gz" | tar -zx && \
    cd ngmlr-* && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make && \
    mv ../bin/ngmlr-*/ngmlr /usr/local/bin/ngmlr && \
    cd ../.. && \
    rm -rf ngmlr-* && \

    # Install seq-align
    git clone --recursive https://github.com/noporpoise/seq-align.git && \
    cd seq-align && \
    make && \
    mv bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf seq-align && \

    # Install STAR
    wget -qO- "https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz" | tar -zx && \
    mv STAR-*/bin/Linux_*_static/* /usr/local/bin/ && \
    rm -rf STAR-* && \

    # Install Unimap
    wget "https://github.com/lh3/unimap/archive/refs/heads/master.zip" && \
    unzip master.zip && \
    cd unimap-master && \
    make && \
    mv unimap /usr/local/bin/unimap && \
    cd .. && \
    rm -rf master.zip unimap-master && \

    # Install wfmash
    wget -qO- "https://github.com/waveygang/wfmash/releases/download/v0.13.0/wfmash-v0.13.0.tar.gz" | tar -zx && \
    cd wfmash-* && \
    cmake -H. -Bbuild -DCMAKE_C_COMPILER="$(which gcc-10)" -DCMAKE_CXX_COMPILER="$(which g++-10)" && \
    cmake --build build -- && \
    mv build/bin/wfmash /usr/local/bin/wfmash && \
    cd .. && \
    rm -rf wfmash-* && \

    # Install Windowmap
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
