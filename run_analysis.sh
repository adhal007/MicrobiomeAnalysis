#!/bin/bash

# Create log directory
mkdir -p logs

# Install system dependencies
echo "Installing system dependencies..."
brew install libxml2
brew install gmp
brew install glpk
brew install libgit2
brew install libssh2
brew install openssl@3

# Install required packages
echo "Installing required packages..."

# First install BiocManager
Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org/')"

# Install igraph first
Rscript -e "install.packages('igraph', repos='https://cloud.r-project.org/')"

# Install Bioconductor packages
Rscript -e "BiocManager::install('dada2', force = TRUE)"
Rscript -e "BiocManager::install('phyloseq', force = TRUE)"

# Install remaining CRAN packages
Rscript -e "install.packages(c('vegan', 'ggplot2', 'pheatmap', 'RColorBrewer', 'ggrepel', 'tidyverse', 'reshape2', 'data.table'), repos='https://cloud.r-project.org/')"

# Install SpiecEasi from GitHub since it's not in Bioconductor 3.21
Rscript -e "if (!requireNamespace('devtools', quietly = TRUE)) install.packages('devtools', repos='https://cloud.r-project.org/')"
Rscript -e "devtools::install_github('zdk123/SpiecEasi')"

# Run the analysis
echo "Running analysis..."
Rscript src/main_analysis.R 2>&1 | tee logs/analysis.log

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo "Analysis completed successfully!"
    echo "Results are available in the results/ directory"
else
    echo "Analysis failed. Please check logs/analysis.log for details"
    exit 1
fi 