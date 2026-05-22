# start from a clean Ubuntu base
FROM ubuntu:22.04

# avoid interactive prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    g++ \
    rsync \
    && rm -rf /var/lib/apt/lists/*

# install Miniconda
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -b -p /opt/conda \
    && rm /tmp/miniconda.sh

ENV PATH /opt/conda/bin:$PATH

# install fastp via conda
RUN conda install -y -c bioconda -c conda-forge fastp \
    && conda clean -afy

# install kraken2 from source
RUN git clone https://github.com/DerrickWood/kraken2.git /tmp/kraken2 \
    && /tmp/kraken2/install_kraken2.sh /usr/local/bin \
    && rm -rf /tmp/kraken2

# install shovill via conda
RUN conda install -y -c bioconda -c conda-forge shovill \
    && conda clean -afy

# install quast via conda
RUN conda install -y -c bioconda -c conda-forge quast \
    && conda clean -afySS

# set working directory inside the container
WORKDIR /data