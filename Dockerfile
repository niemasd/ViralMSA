# Minimal Docker image for ViralMSA using Debian base
FROM debian:stable-slim
RUN apt-get update && \
    apt-get install -y g++ gcc libboost-date-time-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev libboost-system-dev libboost-thread-dev make python3 unzip wget zlib1g-dev && \

    # install ViralMSA v1.1.49
    wget -qO /usr/local/bin/ViralMSA.py "https://github.com/niemasd/ViralMSA/releases/download/1.1.49/ViralMSA.py" && \

    # install Bowtie2 v2.5.5
    wget -q "https://github.com/BenLangmead/bowtie2/releases/download/v2.5.5/bowtie2-2.5.5-linux-x86_64.zip" && \
    unzip bowtie2-*.zip && \
    mv bowtie2-*/bowtie2* /usr/local/bin/ && \
    rm -rf bowtie2-* && \

    # install BWA v0.7.19
    wget -qO- "https://github.com/lh3/bwa/archive/refs/tags/v0.7.19.tar.gz" | tar -zx && \
    cd bwa-* && \
    make && \
    mv bwa bwakit qualfa2fq.pl xa2multi.pl /usr/local/bin/ && \
    cd .. && \
    rm -rf bwa-* && \

    # install DRAGMAP v1.3.0
    wget -qO- "https://github.com/Illumina/DRAGMAP/archive/refs/tags/1.3.0.tar.gz" | tar -zx && \
    cd DRAGMAP-* && \
    HAS_GTEST=0 CFLAGS="-include stdint.h -Wno-error=format-overflow -Wno-error=nonnull" CXXFLAGS="-include cstdint -Wno-error=nonnull" make install && \
    cd .. && \
    rm -rf DRAGMAP-* && \

    # install HISAT2 v2.2.3
    wget -qO- "https://github.com/DaehwanKimLab/hisat2/archive/refs/tags/v2.2.3.tar.gz" | tar -zx && \
    cd hisat2-* && \
    make && \
    mv hisat2 hisat2-* /usr/local/bin/ && \
    cd .. && \
    rm -rf hisat2-*
    # TODO continue here: https://github.com/niemasd/ViralMSA/blob/599bcd8e76023ad71b5da49fac6817d6b24d66df/ViralMSA.py#L794
