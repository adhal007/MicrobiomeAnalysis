# Load required packages
library(tidyverse)
library(data.table)

# Function to check and create directories
create_directories <- function() {
    dirs <- c("data/raw",
              "data/processed",
              "results/figures",
              "results/tables")
    
    for (dir in dirs) {
        if (!dir.exists(dir)) {
            dir.create(dir, recursive = TRUE)
        }
    }
}

# Function to validate input files
validate_input_files <- function(path) {
    # Check if path exists
    if (!dir.exists(path)) {
        stop("Input directory does not exist")
    }
    
    # Check for FASTQ files
    fastq_files <- list.files(path, pattern = "\\.fastq\\.gz$")
    if (length(fastq_files) == 0) {
        stop("No FASTQ files found in input directory")
    }
    
    # Check for metadata file
    if (!file.exists("data/metadata.csv")) {
        stop("Metadata file not found")
    }
    
    return(TRUE)
}

# Function to prepare metadata
prepare_metadata <- function(metadata_path) {
    # Read metadata
    metadata <- read.csv(metadata_path, row.names = 1)
    
    # Add group information based on sample names
    metadata$Group <- ifelse(grepl("a", rownames(metadata)), "Group A", "Group B")
    
    return(metadata)
}

# Function to save analysis results
save_results <- function(results, filename) {
    # Create results directory if it doesn't exist
    dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
    
    # Save results
    write.csv(results, file.path("results/tables", filename))
}

# Function to generate report
generate_report <- function() {
    # Create report directory
    dir.create("results/report", recursive = TRUE, showWarnings = FALSE)
    
    # Generate markdown report
    report <- c(
        "# Microbiome Analysis Report",
        "",
        "## Analysis Summary",
        "",
        "### Alpha Diversity",
        "Alpha diversity metrics were calculated for each sample group.",
        "",
        "### Beta Diversity",
        "Beta diversity was analyzed using PCoA and PERMANOVA.",
        "",
        "### Taxonomic Composition",
        "Relative abundance was analyzed at phylum, genus, and species levels.",
        "",
        "### Network Analysis",
        "Microbial networks were constructed and clusters were identified.",
        "",
        "## Results",
        "",
        "Results are available in the following directories:",
        "- `results/figures/`: Contains all generated plots",
        "- `results/tables/`: Contains statistical results and data tables"
    )
    
    # Write report
    writeLines(report, "results/report/analysis_report.md")
}

# Function to check package dependencies
check_packages <- function() {
    required_packages <- c(
        "dada2", "phyloseq", "vegan", "DESeq2",
        "ggplot2", "pheatmap", "RColorBrewer",
        "ggpubr", "ggrepel", "igraph", "SpiecEasi",
        "tidyverse", "reshape2", "data.table"
    )
    
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
    
    if (length(missing_packages) > 0) {
        stop(paste("The following packages are missing:", paste(missing_packages, collapse = ", ")))
    }
    
    return(TRUE)
}

# Function to clean up temporary files
cleanup_temp_files <- function() {
    temp_dirs <- c("data/raw/filtered")
    
    for (dir in temp_dirs) {
        if (dir.exists(dir)) {
            unlink(dir, recursive = TRUE)
        }
    }
} 