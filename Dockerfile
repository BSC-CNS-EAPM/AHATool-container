################## BASE IMAGE ######################

FROM ubuntu:latest

################## METADATA ######################

LABEL maintainer="Albert Cañellas-Sole <albert.canellas@bsc.es>" \
    container="AHATool" \
    about.summary="AHATool an Automatic HMM and Analysis Tool" \
    about.home="https://github.com/BSC-CNS-EAPM/AHATool-container" \
    software.version="1.0"

################## INSTALLATION ######################

# Update to latest packages
RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git zip libz-dev byobu samtools libncursesw5-dev python3-pip libbz2-dev lzma-dev liblzma-dev \
    libcurl4-gnutls-dev ncbi-blast+ perl libxml-simple-perl cpanminus libwww-perl libnet-perl \
    libssl-dev libio-socket-ssl-perl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add AHATool resources
WORKDIR /home/AHATool/AHATool_Resources
RUN wget https://github.com/BSC-CNS-EAPM/AHATool-container/raw/main/AHATool_Resources/shflags && \
    wget https://github.com/BSC-CNS-EAPM/AHATool-container/raw/main/AHATool_Resources/update_FASTAdb.pl && \
    wget https://github.com/BSC-CNS-EAPM/AHATool-container/raw/main/AHATool_Resources/SOFTWAREneeded.txt 

# Install Signalp6
ADD ./AHATool_Resources/signalp6.tar .
WORKDIR /home/AHATool/AHATool_Resources/signalp6_fast
RUN pip install signalp-6-package/ 
RUN cp -r signalp-6-package/models/* $(python3 -c "import signalp; import os; print(os.path.dirname(signalp.__file__))")/model_weights/

# Install EDirect
WORKDIR /home/AHATool/AHATool_Resources
RUN sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
ENV PATH=$PATH:/root/edirect

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
    cd easel && make install && \
    rm /home/AHATool/AHATool_Resources/hmmer/hmmer.tar.gz

# Install T_coffe
WORKDIR /home/AHATool/AHATool_Resources/t_coffee
RUN wget https://s3.eu-central-1.amazonaws.com/tcoffee-packages/Stable/Latest/T-COFFEE_distribution_Version_13.45.0.4846264.tar.gz && \
    tar xvf T-COFFEE_distribution_Version_13.45.0.4846264.tar.gz && \
    cd T-COFFEE_distribution_Version_13.45.0.4846264 && \ 
    ./install tcoffee && \
    rm /home/AHATool/AHATool_Resources/t_coffee/T-COFFEE_distribution_Version_13.45.0.4846264.tar.gz
ENV PATH=$PATH:/root/.t_coffee/bin/linux

# Add AHATool
WORKDIR /home/AHATool
RUN wget https://github.com/BSC-CNS-EAPM/AHATool-container/raw/main/AHATool.sh && \
    chmod +x AHATool.sh
ENV PATH=$PATH:/home/AHATool

WORKDIR /home/projects
 






