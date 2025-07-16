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
        ca-certificates \
        r-base \
        r-base-dev \
        libxml2-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        git \
        build-essential \
        gfortran \
        libblas-dev \
        liblapack-dev \
        libgsl-dev && \
    apt-get clean

# download liftOver binary and make executable
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod +x liftOver && \
    mv liftOver /usr/local/bin/

# install R dependencies & ichorCNA from github
# clone ichorCNA repo
RUN git clone https://github.com/broadinstitute/ichorCNA.git /opt/ichorCNA
# make ichorCNA scripts available globally
ENV PATH="/opt/ichorCNA/scripts:$PATH"
# install R dependencies step by step with better error handling
RUN Rscript -e "install.packages(c('optparse', 'data.table', 'foreach', 'doParallel', 'ggplot2', 'gridExtra', 'remotes'), repos='https://cloud.r-project.org')"
# install BiocManager
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"
# install Bioconductor dependencies one by one to isolate any issues
RUN Rscript -e "BiocManager::install('S4Vectors', update=FALSE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('IRanges', update=FALSE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('XVector', update=FALSE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('GenomeInfoDb', update=FALSE, ask=FALSE)"
# try installing GenomicRanges with increased memory and timeout
RUN Rscript -e "options(timeout=300); BiocManager::install('GenomicRanges', update=FALSE, ask=FALSE, force=TRUE)"
# install HMMcopy from Bioconductor (more reliable than CRAN version)
RUN Rscript -e "BiocManager::install('HMMcopy', update=FALSE, ask=FALSE)"
# install ichorCNA from local repo
RUN Rscript -e "remotes::install_local('/opt/ichorCNA', dependencies=TRUE, upgrade='never')"
# verify installation
RUN Rscript -e "library(ichorCNA); cat('ichorCNA successfully installed\\n')"
# debug: print all installed packages
RUN Rscript -e "installed <- rownames(installed.packages()); writeLines(installed)"