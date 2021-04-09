# ViralMSA minimal Docker image using
FROM python:3.8-slim-buster
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# Set up environment and install dependencies
RUN apt-get update && apt-get -y upgrade && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y bzip2 git gsl-bin libgsl0-dev perl unzip wget zlib1g-dev && \
    pip3 install biopython

# Install bowtie2 (v.2.4.1)
RUN wget "https://github.com/BenLangmead/bowtie2/releases/download/v2.4.1/bowtie2-2.4.1-linux-x86_64.zip" && \
    unzip bowtie2-*.zip && mv bowtie2-2.4.1-linux-x86_64/bowtie2* /usr/local/bin && rm -rf bowtie2-*

# Install HISAT2 (2.2.1)
RUN wget -O hisat2.zip "https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download" && \
    unzip hisat2.zip && mv hisat2-*/hisat2* /usr/local/bin && rm -rf hisat2*

# Install Minimap2 (2.17)
RUN wget -qO- "https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2" | tar xj && \
    mv minimap2-*/minimap2 /usr/local/bin && rm -rf minimap2-*

# Install STAR (2.7.5c)
RUN wget -qO- "https://github.com/alexdobin/STAR/archive/2.7.5c.tar.gz" | tar -zx && \
    mv STAR-*/bin/Linux_x86_64_static/* /usr/local/bin && rm -rf STAR-*

# Install wfmash
#RUN git clone https://github.com/ekg/wfmash.git && \
#    cd wfmash && ./bootstrap.sh && ./configure && make && make install && cd .. && rm -rf wfmash

# Set up ViralMSA
RUN wget -O /usr/local/bin/ViralMSA.py "https://raw.githubusercontent.com/niemasd/ViralMSA/master/ViralMSA.py" && chmod a+x /usr/local/bin/ViralMSA.py

# Clean up
RUN rm -rf /root/.cache && \
    rm -rf /tmp/*
