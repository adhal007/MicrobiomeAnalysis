# Microbiome Analysis Pipeline

This repository contains a comprehensive pipeline for analyzing 16S rRNA sequencing data using DADA2 and various statistical methods for microbiome analysis.

## Project Structure

```
.
├── data/                  # Raw and processed data
├── src/                   # Source code
│   ├── main_analysis.R    # Main analysis script
│   ├── visualization.R    # Visualization functions
│   └── utils.R           # Utility functions
├── results/              # Analysis outputs
│   ├── figures/          # Generated plots
│   └── tables/           # Statistical results
├── requirements.txt      # R package dependencies
└── README.md            # This file
```

## Analysis Pipeline

The analysis pipeline includes:

1. Quality control and preprocessing using DADA2
2. Taxonomic assignment
3. Alpha diversity analysis
   - Shannon index
   - Simpson index
   - Chao1 index
4. Beta diversity analysis
   - PCoA plots
   - PERMANOVA tests
5. Compositional analysis
   - Heatmaps at phylum, genus, and species levels
6. Microbial network analysis
   - Co-occurrence networks
   - Cluster identification

## Statistical Analysis

The pipeline performs statistical comparisons between groups:
- Group A (samples with 'a' in name) vs Group B (samples with 'b' in name)
- Statistical tests include:
  - PERMANOVA for beta diversity
  - Kruskal-Wallis test for alpha diversity
  - Differential abundance analysis

## Requirements

See `requirements.txt` for the complete list of R packages required.

## Usage

1. Place your FASTQ files in the `data/raw` directory
2. Update the sample metadata in `data/metadata.csv`
3. Run the analysis:
   ```R
   source("src/main_analysis.R")
   ```

## Output

The analysis generates:
- Taxonomic tables at different levels
- Diversity metrics
- Statistical test results
- Publication-ready figures
- Network analysis results 