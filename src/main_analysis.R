# Load required packages
library(dada2)
library(phyloseq)
library(vegan)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(igraph)
library(SpiecEasi)
library(tidyverse)
library(reshape2)
library(data.table)

# Set working directory and create output directories
setwd("..")
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# Function to process FASTQ files with DADA2
process_fastq_files <- function(path) {
    # Forward and reverse fastq filenames
    fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
    fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
    
    # Extract sample names
    sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
    
    # Quality filtering
    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
    
    # Filter and trim
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                        truncLen=c(210,160),
                        maxN=0, maxEE=c(2,2), truncQ=2,
                        rm.phix=TRUE,
                        compress=TRUE, multithread=TRUE)
    
    # Learn error rates
    errF <- learnErrors(filtFs, multithread=TRUE)
    errR <- learnErrors(filtRs, multithread=TRUE)
    
    # Dereplicate reads
    derepFs <- derepFastq(filtFs, verbose=TRUE)
    derepRs <- derepFastq(filtRs, verbose=TRUE)
    
    # Sample inference
    dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
    dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
    
    # Merge paired reads
    mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
    
    # Construct sequence table
    seqtab <- makeSequenceTable(mergers)
    
    # Remove chimeras
    seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                       multithread=TRUE, verbose=TRUE)
    
    # Track reads through pipeline
    getN <- function(x) sum(getUniques(x))
    track <- cbind(out, 
                  sapply(dadaFs, getN),
                  sapply(dadaRs, getN),
                  sapply(mergers, getN),
                  rowSums(seqtab.nochim))
    colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
    rownames(track) <- sample.names
    
    # Save tracking information
    write.csv(track, "results/tables/read_tracking.csv")
    
    return(seqtab.nochim)
}

# Function to perform taxonomic assignment
assign_taxonomy <- function(seqtab) {
    # Assign taxonomy using Silva database
    taxa <- assignTaxonomy(seqtab, 
                          "data/tax/silva_nr_v132_train_set.fa.gz",
                          multithread = TRUE)
    
    # Add species level assignment
    taxa <- addSpecies(taxa, "data/tax/silva_species_assignment_v132.fa.gz")
    
    # Save taxonomic assignments
    write.csv(taxa, "results/tables/taxonomic_assignments.csv")
    
    return(taxa)
}

# Function to create metadata from sample names
create_metadata <- function(sample_names) {
    metadata <- data.frame(
        Sample = sample_names,
        Group = ifelse(grepl("a", sample_names), "Group A", "Group B"),
        row.names = sample_names
    )
    write.csv(metadata, "results/tables/sample_metadata.csv")
    return(metadata)
}

# Function to filter out contaminants and unwanted taxa
filter_contaminants <- function(taxa) {
    # Convert taxa matrix to data frame
    tax_table <- as.data.frame(taxa)
    
    # Print initial statistics
    cat("Initial status:\n")
    cat("Total taxa:", nrow(tax_table), "\n")
    
    # Print kingdom-level distribution
    kingdom_counts <- table(tax_table$Kingdom)
    cat("\nKingdom distribution:\n")
    print(kingdom_counts)
    
    # Keep only bacteria (remove eukaryotes, chloroplasts, mitochondria)
    keep_taxa <- rownames(tax_table)[
        (!is.na(tax_table$Kingdom) & tax_table$Kingdom == "Bacteria") &
        !grepl("Chloroplast|Mitochondria", tax_table$Family, ignore.case = TRUE) &
        !grepl("Chloroplast|Mitochondria", tax_table$Order, ignore.case = TRUE)
    ]
    
    # Check if we have any taxa left
    if(length(keep_taxa) == 0) {
        stop("No bacterial taxa found after filtering!")
    }
    
    # Print number of taxa retained
    cat("\nTaxa after filtering:", length(keep_taxa), "\n")
    
    # Print family-level distribution for retained taxa
    family_counts <- table(tax_table[keep_taxa, "Family"])
    cat("\nFamily distribution of retained taxa:\n")
    print(family_counts)
    
    # Create filtered taxonomy table
    taxa_filtered <- taxa[keep_taxa, ]
    
    # Save filtering statistics
    filter_stats <- data.frame(
        Metric = c("Original taxa", "Filtered taxa",
                  "Number of families"),
        Value = c(nrow(taxa), nrow(taxa_filtered),
                 length(unique(na.omit(tax_table[keep_taxa, "Family"]))))
    )
    write.csv(filter_stats, "results/tables/taxonomy_filtering_statistics.csv")
    
    return(taxa_filtered)
}

# Function to create phyloseq object
create_phyloseq <- function(seqtab, taxa, sample_names) {
    # Create metadata
    metadata <- create_metadata(sample_names)
    
    # Filter taxonomy to keep only bacterial families
    taxa_filtered <- filter_contaminants(taxa)
    
    # Subset sequence table to keep only filtered taxa
    seqtab_filtered <- seqtab[, colnames(seqtab) %in% rownames(taxa_filtered)]
    
    # Create phyloseq object
    ps <- phyloseq(otu_table(seqtab_filtered, taxa_are_rows = FALSE),
                   tax_table(taxa_filtered),
                   sample_data(metadata))
    
    # Remove any samples with zero reads after filtering
    ps <- prune_samples(sample_sums(ps) > 0, ps)
    
    # Print final statistics
    cat("\nFinal phyloseq object status:\n")
    cat("Samples:", nsamples(ps), "\n")
    cat("Taxa:", ntaxa(ps), "\n")
    cat("Total reads:", sum(sample_sums(ps)), "\n")
    
    return(ps)
}

# Function to calculate alpha diversity
calculate_alpha_diversity <- function(ps) {
    # Calculate various alpha diversity metrics
    alpha_div <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Chao1"))
    
    # Add sample information
    alpha_div$Sample <- rownames(alpha_div)
    alpha_div$Group <- sample_data(ps)$Group
    
    # Perform statistical tests
    shannon_test <- wilcox.test(Shannon ~ Group, data = alpha_div)
    simpson_test <- wilcox.test(Simpson ~ Group, data = alpha_div)
    chao1_test <- wilcox.test(Chao1 ~ Group, data = alpha_div)
    
    # Save test results
    test_results <- data.frame(
        Metric = c("Shannon", "Simpson", "Chao1"),
        p_value = c(shannon_test$p.value, simpson_test$p.value, chao1_test$p.value)
    )
    write.csv(test_results, "results/tables/alpha_diversity_tests.csv")
    
    return(alpha_div)
}

# Function to calculate beta diversity
calculate_beta_diversity <- function(ps) {
    # Remove samples with zero reads
    ps <- prune_samples(sample_sums(ps) > 0, ps)
    
    # Remove taxa with zero reads
    ps <- prune_taxa(taxa_sums(ps) > 0, ps)
    
    # Convert to relative abundance to handle different sequencing depths
    ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
    
    # Replace any remaining zeros with a small value to avoid infinite distances
    otu_table(ps_rel)[otu_table(ps_rel) == 0] <- 1e-10
    
    # Calculate distance matrix
    dist_matrix <- distance(ps_rel, method = "bray")
    
    # Check for any infinite or NA values
    if(any(is.infinite(dist_matrix)) || any(is.na(dist_matrix))) {
        warning("Infinite or NA values detected in distance matrix. Using alternative method.")
        # Use alternative method: Hellinger transformation
        ps_hell <- transform_sample_counts(ps, function(x) sqrt(x/sum(x)))
        dist_matrix <- distance(ps_hell, method = "bray")
    }
    
    # Perform PCoA
    pcoa <- ordinate(ps_rel, method = "PCoA", distance = "bray")
    
    # Prepare metadata for PERMANOVA
    metadata <- data.frame(Group = sample_data(ps)$Group)
    rownames(metadata) <- sample_names(ps)
    
    # Perform PERMANOVA
    permanova <- adonis2(dist_matrix ~ Group, data = metadata)
    write.csv(permanova, "results/tables/permanova_results.csv")
    
    # Calculate and save dispersion test results
    disp <- betadisper(dist_matrix, metadata$Group)
    disp_test <- permutest(disp)
    capture.output(disp_test, file = "results/tables/beta_dispersion_test.txt")
    
    # Save distance matrix
    write.csv(as.matrix(dist_matrix), "results/tables/beta_diversity.csv")
    
    return(list(dist_matrix = dist_matrix, pcoa = pcoa))
}

# Function to create composition heatmaps
create_composition_heatmap <- function(ps, level) {
    # Aggregate to specified taxonomic level
    ps_agg <- tax_glom(ps, taxrank = level)
    
    # Convert to relative abundance
    ps_agg_rel <- transform_sample_counts(ps_agg, function(x) x/sum(x))
    
    # Extract abundance matrix
    abund_matrix <- as.matrix(otu_table(ps_agg_rel))
    
    # Get top taxa
    top_taxa <- names(sort(rowSums(abund_matrix), decreasing = TRUE)[1:20])
    abund_matrix <- abund_matrix[top_taxa, ]
    
    # Create heatmap
    pheatmap(abund_matrix,
             scale = "row",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             show_rownames = TRUE,
             show_colnames = TRUE,
             filename = paste0("results/figures/heatmap_", level, ".pdf"))
}

# Function to perform network analysis
perform_network_analysis <- function(ps) {
    # Convert to relative abundance
    ps_rel <- transform_sample_counts(ps, function(x) x/sum(x))
    
    # Extract abundance matrix
    abund_matrix <- as.matrix(otu_table(ps_rel))
    
    # Remove taxa with zero variance
    var_taxa <- apply(abund_matrix, 1, var)
    abund_matrix <- abund_matrix[var_taxa > 0, ]
    
    # Check if we have any taxa left
    if (nrow(abund_matrix) == 0) {
        warning("No taxa with non-zero variance found. Skipping network analysis.")
        return(list(network = NULL, clusters = NULL))
    }
    
    # Check if we have enough taxa for correlation
    if (nrow(abund_matrix) < 2) {
        warning("Insufficient number of taxa for network analysis. Need at least 2 taxa.")
        return(list(network = NULL, clusters = NULL))
    }
    
    # Calculate correlation matrix with handling of NaN values
    cor_matrix <- cor(t(abund_matrix), method = "spearman", use = "pairwise.complete.obs")
    
    # Replace any remaining NaN values with 0
    cor_matrix[is.na(cor_matrix)] <- 0
    
    # Create network with absolute correlations
    net <- graph.adjacency(abs(cor_matrix), mode = "undirected", weighted = TRUE, diag = FALSE)
    
    # Remove edges with zero weight
    net <- delete.edges(net, which(E(net)$weight == 0))
    
    # Check if we have any edges left
    if (ecount(net) == 0) {
        warning("No significant correlations found. Network is empty.")
        return(list(network = NULL, clusters = NULL))
    }
    
    # Identify clusters
    clusters <- cluster_fast_greedy(net)
    
    # Save network metrics
    network_metrics <- data.frame(
        Nodes = vcount(net),
        Edges = ecount(net),
        Density = graph.density(net),
        Modularity = modularity(clusters)
    )
    write.csv(network_metrics, "results/tables/network_metrics.csv")
    
    return(list(network = net, clusters = clusters))
}

# Function to plot beta diversity results
plot_beta_diversity <- function(beta_div, ps) {
    # Extract PCoA coordinates from the vectors component
    pcoa_scores <- beta_div$pcoa$vectors
    
    # Get sample names from PCoA scores
    pcoa_samples <- rownames(pcoa_scores)
    
    # Subset phyloseq object to match PCoA samples
    ps_subset <- prune_samples(pcoa_samples, ps)
    
    # Create data frame for plotting
    plot_data <- data.frame(
        PC1 = pcoa_scores[, 1],
        PC2 = pcoa_scores[, 2],
        Group = sample_data(ps_subset)$Group,
        Sample = pcoa_samples
    )
    
    # Calculate variance explained
    var_explained <- round(beta_div$pcoa$values$Relative_eig * 100, 1)
    
    # Create PCoA plot
    p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3) +
        stat_ellipse(level = 0.95, type = "t") +  # Use t-distribution for ellipses
        geom_text_repel(aes(label = Sample), 
                       size = 3, 
                       max.overlaps = 50,  # Increase max overlaps
                       min.segment.length = 0.1,  # Reduce minimum segment length
                       box.padding = 0.5) +  # Add padding around labels
        scale_color_brewer(palette = "Set1") +
        labs(
            x = paste0("PC1 (", var_explained[1], "%)"),
            y = paste0("PC2 (", var_explained[2], "%)"),
            title = "Principal Coordinates Analysis (PCoA)",
            subtitle = "Based on Bray-Curtis dissimilarity"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right"
        )
    
    # Save the plot
    ggsave("results/figures/beta_diversity_pcoa.pdf", p, width = 12, height = 10)  # Increase plot size
    
    # Create distance boxplot
    dist_matrix <- as.matrix(beta_div$dist_matrix)
    
    # Use reshape2::melt explicitly to avoid conflicts
    dist_data <- reshape2::melt(dist_matrix)
    colnames(dist_data) <- c("Sample1", "Sample2", "Distance")
    
    # Remove self-comparisons and duplicates
    dist_data <- dist_data[dist_data$Sample1 < dist_data$Sample2, ]
    
    # Add group information
    dist_data$Group1 <- sample_data(ps_subset)[dist_data$Sample1, "Group"]
    dist_data$Group2 <- sample_data(ps_subset)[dist_data$Sample2, "Group"]
    dist_data$Comparison <- ifelse(dist_data$Group1 == dist_data$Group2, 
                                 "Within Group", "Between Groups")
    
    # Create boxplot
    p_box <- ggplot(dist_data, aes(x = Comparison, y = Distance, fill = Comparison)) +
        geom_boxplot() +
        scale_fill_brewer(palette = "Set2") +
        labs(
            title = "Beta Diversity Distances",
            subtitle = "Bray-Curtis dissimilarity",
            x = "Comparison Type",
            y = "Distance"
        ) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "none"
        )
    
    # Save the boxplot
    ggsave("results/figures/beta_diversity_boxplot.pdf", p_box, width = 8, height = 6)
    
    return(list(pcoa_plot = p, boxplot = p_box))
}

# Main analysis pipeline
main <- function() {
    # Process FASTQ files
    seqtab <- process_fastq_files("data/raw")
    
    # Get sample names from sequence table
    sample_names <- rownames(seqtab)
    
    # Assign taxonomy
    taxa <- assign_taxonomy(seqtab)
    
    # Create phyloseq object
    ps <- create_phyloseq(seqtab, taxa, sample_names)
    
    # Calculate alpha diversity
    alpha_div <- calculate_alpha_diversity(ps)
    
    # Save alpha diversity results
    write.csv(alpha_div, "results/tables/alpha_diversity.csv")
    
    # Calculate beta diversity
    beta_div <- calculate_beta_diversity(ps)
    
    # Save beta diversity results
    write.csv(as.matrix(beta_div$dist_matrix), "results/tables/beta_diversity.csv")
    
    # Plot beta diversity results
    beta_plots <- plot_beta_diversity(beta_div, ps)
    
    # Create composition heatmaps
    for (level in c("Phylum", "Genus", "Species")) {
        create_composition_heatmap(ps, level)
    }
    
    # Perform network analysis
    network_results <- perform_network_analysis(ps)
    
    # Save network results
    saveRDS(network_results, "results/tables/network_analysis.rds")
}

# # Run the analysis
# main() 
