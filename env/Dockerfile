# =====================================================================
#  josephhalstead/dragenwgs_post_processing Dockerfile
#
# =====================================================================

# Use miniconda parent image - debian based
FROM continuumio/miniconda:4.7.10
LABEL maintainer="Joseph Halstead"

# Set the working directory & copy across current directory
# Set python so that it doesn't create .pyc files
WORKDIR /var/app
COPY . /var/app
ENV PYTHONDONTWRITEBYTECODE=true

# Install awk and ps for use with Nextflow reports
# remove apt cache after running to reduce image size
RUN apt-get update \
 && apt-get install --yes --no-install-recommends \
        build-essential \
        gawk \
        git \
        procps \
        zlib1g-dev \
        libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

# Update conda, make the conda environment, clean unnecessary files
# afterwards
RUN conda install \
        python=3.7 \
        bioconda::bcftools=1.9 \
        bioconda::vt=2015.11.10 \
        bioconda::vcftools=0.1.16 \
        bioconda::ensembl-vep=100 \
        bioconda::nextflow=19 \
        bioconda::slivar=0.1.11 \
        bioconda::pysam=0.15.2

RUN pip install pandas==0.23.4
RUN pip install pyvariantfilter==1.0.4
RUN pip install pyupdog==1.0.2

