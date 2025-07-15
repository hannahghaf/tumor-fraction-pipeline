# ubuntu 22.04 as base image
# build image for x86_64 architecture for liftOver
FROM --platform=linux/amd64 ubuntu:22.04

# prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# update apt package list, install tools
RUN apt-get update && \
    apt-get install -y \
        bedtools \
        samtools \
        wget \
        unzip \
        libcurl4 \
        ca-certificates && \
    apt-get clean

# download liftOver binary and make executable
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod +x liftOver && \
    mv liftOver /usr/local/bin/