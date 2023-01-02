################## BASE IMAGE ######################

FROM ubuntu:latest

################## METADATA ######################

LABEL maintainer="Albert Cañellas-Sole <albert.canellas@bsc.es>" \
    container="AHATool" \
    about.summary="design catalytic residues, perform in silico directed evolution of an existing active site." \
    about.home="https://github.com/BSC-CNS-EAPM/AHATool-container" \
    software.version="1.0"

################## INSTALLATION ######################

# Update to latest packages
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git zip libz-dev byobu samtools libncursesw5-dev libbz2-dev lzma-dev liblzma-dev \
    libcurl4-gnutls-dev hmmer ncbi-blast+ perl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /home/AHATool/AHATool_Resources
ADD ./AHATool_Resources/shflags .
ADD ./AHATool_Resources/update_FASTAdb.pl .

WORKDIR /home/AHATool
ADD AHATool.sh .

 






