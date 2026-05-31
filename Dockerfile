# Dockerfile
FROM rocker/r-ver:4.3.0

# Install Linux system dependencies 
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgsl-dev \
    libglpk-dev \
    libpng-dev \
    zlib1g-dev \
    make \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Fix CRAN repository configuration and install BiocManager
RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages('BiocManager')"

# Install explicit workflow packages
RUN R -e "BiocManager::install(c('Seurat', 'gprofiler2', 'ggplot2', 'dplyr', 'magrittr', 'Matrix', 'limma'))"

# Initialize isolated execution workspace
WORKDIR /workspace

# Run complete orchestrated pipeline by default
CMD ["Rscript", "main.R"]