FROM ubuntu:latest

# Update and install basic dependencies
# findbin a hidden dep for barrnap
RUN apt-get update && apt-get install -y wget ca-certificates

# Make a directory to hold all our new tools and an updated PATH variable to use them
RUN mkdir -p /opt/bin
ENV PATH="/opt/bin:${PATH}"
WORKDIR /opt

# Download bedtools and make it executable
RUN wget -q -O ./bin/bedtools https://github.com/arq5x/bedtools2/releases/download/v2.29.2/bedtools.static.binary && \
        chmod ugo+x ./bin/bedtools

# Download and unpack BLAST
RUN wget -q https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.1/ncbi-blast-2.10.1+-x64-linux.tar.gz && \
        tar -zxvf ncbi-blast-2.10.1+-x64-linux.tar.gz && \
        ln -s /opt/ncbi-blast-2.10.1+/bin/* /opt/bin

# Download and unpack barrnap
RUN wget -q -O barrnap-0.9.tar.gz https://github.com/tseemann/barrnap/archive/0.9.tar.gz && \
        tar -zxvf barrnap-0.9.tar.gz && \
        ln -s /opt/barrnap-0.9/bin/barrnap /opt/bin

# Download and unpack salmon
RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.3.0/salmon-1.3.0_linux_x86_64.tar.gz && \
        tar -zxvf salmon-1.3.0_linux_x86_64.tar.gz && \
        ln -s /opt/salmon-latest_linux_x86_64/bin/salmon /opt/bin

# Install bbtools
RUN apt-get install -y openjdk-8-jre && \
        wget -q -O bbtools.tar.gz https://sourceforge.net/projects/bbmap/files/latest/download && \
        tar -zxvf bbtools.tar.gz && \
        ln -s /opt/bbmap/bbduk.sh /opt/bin/bbduk.sh && \
        ln -s /opt/bbmap/dedupe.sh /opt/bin/dedupe.sh

# We get issues if we don't install this after jre (hidden dep needed by barrnap)
RUN apt-get install -y libfindbin-libs-perl

# Install gffread
RUN wget -q https://github.com/gpertea/gffread/releases/download/v0.12.3/gffread-0.12.3.Linux_x86_64.tar.gz && \
        tar -zxvf gffread-0.12.3.Linux_x86_64.tar.gz && \
        cp gffread-0.12.3.Linux_x86_64/gffread /opt/bin

RUN useradd -ms /bin/bash ribocleaner

USER ribocleaner

RUN echo "umask 002" >> /home/ribocleaner/.bashrc
