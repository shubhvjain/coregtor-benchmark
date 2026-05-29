FROM rocker/r-ver:4.3.0

# Install system dependencies required by TCGAbiolinks
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    make \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager and TCGAbiolinks
RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('TCGAbiolinks', 'SummarizedExperiment'))"

# Set up directories
RUN mkdir /app /data
# COPY download_tcga.R /app/download_tcga.R
# RUN chmod +x /app/download_tcga.R

WORKDIR /app

# Set the script to run automatically on container startup
# ENTRYPOINT ["Rscript", "/app/setup/download_tcga.R"]