# Minimal Docker image for ViralMSA using Debian base
FROM debian:stable-slim
RUN apt-get update && \
    apt-get install -y bzip2 g++ gcc git libboost-date-time-dev libboost-filesystem-dev libboost-iostreams-dev libboost-program-options-dev libboost-system-dev libboost-thread-dev libbz2-dev libcurl4-openssl-dev liblzma-dev make python3 python-is-python3 unzip wget zlib1g-dev && \
    mkdir -p /usr/include/sys && \
    ln -s /usr/include/linux/sysctl.h /usr/include/sys/sysctl.h && \

    # install ViralMSA v1.1.49
    wget -qO /usr/local/bin/ViralMSA.py "https://github.com/niemasd/ViralMSA/releases/download/1.1.49/ViralMSA.py" && \

    # install htslib v1.24
    wget -qO- "https://github.com/samtools/htslib/releases/download/1.24/htslib-1.24.tar.bz2" | tar -xj && \
    cd htslib-* && \
    ./configure && \
    make install && \
    cd .. && \
    rm -rf htslib-* && \
    ldconfig && \

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
    rm -rf hisat2-* && \

    # install LRA v1.3.7.2
    wget -qO- "https://github.com/ChaissonLab/LRA/archive/refs/tags/1.3.7.2.tar.gz" | tar -zx && \
    cd LRA-* && \
    make && \
    mv lra /usr/local/bin/ && \
    cd .. && \
    rm -rf LRA-* && \

    # install Minimap2 v2.31
    wget -qO- "https://github.com/lh3/minimap2/releases/download/v2.31/minimap2-2.31_x64-linux.tar.bz2" | tar -jx && \
    mv minimap2-*/minimap2 /usr/local/bin/ && \
    rm -rf minimap2-* && \

    # install mm2-fast v1.0
    wget -qO- "https://github.com/bwa-mem2/mm2-fast/releases/download/mm2-fast-v1.0/Source_code_including_submodules.tar.gz" | tar -zx && \
    cd mm2-fast && \
    make && \
    mv minimap2 /usr/local/bin/mm2-fast && \
    cd .. && \
    rm -rf mm2-fast && \

    # install NGMLR v0.2.7
    wget -qO- "https://github.com/philres/ngmlr/releases/download/v0.2.7/ngmlr-0.2.7-linux-x86_64.tar.gz" | tar -zx && \
    mv ngmlr-*/ngmlr /usr/local/bin/ && \
    rm -rf ngmlr-* && \

    # install Rammap v1.2.0
    wget -qO /usr/local/bin/rammap "https://github.com/jwanglab/rammap/releases/download/v1.2.0/rammap_x86_64-unknown-linux-gnu_v1.2.0" && \
    chmod a+x /usr/local/bin/rammap && \

    # install seq-align
    git clone --recurse-submodules https://github.com/noporpoise/seq-align.git && \
    cd seq-align && \
    make && \
    mv bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf seq-align && \

    # install STAR v2.7.11b
    wget -q "https://github.com/alexdobin/STAR/releases/download/2.7.11b/STAR_2.7.11b.zip" && \
    unzip STAR_*.zip && \
    mv STAR_*/Linux_x86_64_static/* /usr/local/bin/ && \
    rm -rf STAR_* && \

    # install unimap
    wget -qO- "https://github.com/lh3/unimap/archive/refs/heads/master.tar.gz" | tar -zx && \
    cd unimap-* && \
    make && \
    mv unimap /usr/local/bin/ && \
    cd .. && \
    rm -rf unimap-* && \

    # install wfmash v0.24.2
    wget -qO /usr/local/bin/wfmash "https://github.com/waveygang/wfmash/releases/download/v0.24.2/wfmash-v0.24.2-linux-x86_64" && \
    chmod a+x /usr/local/bin/wfmash && \

    # install Winnowmap v2.03
    wget -qO- "https://github.com/marbl/Winnowmap/archive/refs/tags/v2.03.tar.gz" | tar -zx && \
    cd Winnowmap-* && \
    make && \
    mv bin/* /usr/local/bin/ && \
    cd .. && \
    rm -rf Winnowmap-*
