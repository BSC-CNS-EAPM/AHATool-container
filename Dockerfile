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
    apt-get install -y wget bzip2 ca-certificates curl git zip libz-dev byobu samtools libncursesw5-dev python3-pip libbz2-dev lzma-dev liblzma-dev \
    libcurl4-gnutls-dev ncbi-blast+ perl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add AHATool resources
WORKDIR /home/AHATool/AHATool_Resources
ADD ./AHATool_Resources/shflags .
ADD ./AHATool_Resources/update_FASTAdb.pl .
ADD ./AHATool_Resources/SOFTWAREneeded.txt .

# Install Signalp6
ADD ./AHATool_Resources/signalp6.tar .
WORKDIR /home/AHATool/AHATool_Resources/signalp6_fast
RUN pip install signalp-6-package/ 
RUN cp -r signalp-6-package/models/* $(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))")/model_weights/

# Install EDirect
WORKDIR /home/AHATool/AHATool_Resources/edirect
RUN sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
RUN echo "export PATH=\${PATH}:/root/edirect" >> ${HOME}/.bashrc

# Install HMMER
WORKDIR /home/AHATool/AHATool_Resources/hmmer
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
    tar zxf hmmer.tar.gz && \
    cd hmmer-3.3.2 && \
    mkdir build && \
    ./configure && \
    make && \
    make check && \
    make install && \
    cd easel && make install

# Add AHATool
WORKDIR /home/AHATool
ADD AHATool.sh .

 






