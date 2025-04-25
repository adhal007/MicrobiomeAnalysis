# Load required packages
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(ggrepel)

# Function to plot alpha diversity
plot_alpha_diversity <- function(alpha_div) {
    # Create plots for each diversity metric
    metrics <- c("Shannon", "Simpson", "Chao1")
    plots <- list()
    
    for (metric in metrics) {
        p <- ggplot(alpha_div, aes(x = Group, y = get(metric), fill = Group)) +
            geom_boxplot(alpha = 0.7) +
            geom_jitter(width = 0.2, size = 2) +
            stat_compare_means(method = "wilcox.test", label = "p.signif") +
            labs(title = paste(metric, "Diversity"),
                 y = metric,
                 x = "Group") +
            theme_minimal() +
            theme(legend.position = "none")
        
        plots[[metric]] <- p
    }
    
    # Combine plots
    combined_plot <- ggarrange(plotlist = plots, ncol = 3)
    ggsave("results/figures/alpha_diversity.pdf", combined_plot, width = 12, height = 4)
}

# Function to plot PCoA
plot_pcoa <- function(pcoa, metadata) {
    # Extract coordinates
    coords <- as.data.frame(pcoa$vectors[, 1:2])
    colnames(coords) <- c("PC1", "PC2")
    coords$Group <- ifelse(grepl("a", rownames(coords)), "Group A", "Group B")
    
    # Calculate variance explained
    var_explained <- round(pcoa$values$Relative_eig * 100, 1)
    
    # Create plot
    p <- ggplot(coords, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3) +
        stat_ellipse(type = "t", linetype = 2) +
        labs(x = paste0("PC1 (", var_explained[1], "%)"),
             y = paste0("PC2 (", var_explained[2], "%)"),
             title = "PCoA Plot") +
        theme_minimal() +
        theme(legend.position = "right")
    
    ggsave("results/figures/pcoa_plot.pdf", p, width = 8, height = 6)
}

# Function to plot taxonomic composition
plot_taxonomic_composition <- function(ps, level) {
    # Aggregate to specified level
    ps_agg <- tax_glom(ps, taxrank = level)
    
    # Convert to relative abundance
    ps_agg_rel <- transform_sample_counts(ps_agg, function(x) x/sum(x))
    
    # Melt data for plotting
    melted_data <- psmelt(ps_agg_rel)
    
    # Get top taxa
    top_taxa <- names(sort(taxa_sums(ps_agg_rel), decreasing = TRUE)[1:10])
    
    # Create plot
    p <- ggplot(melted_data, aes(x = Sample, y = Abundance, fill = get(level))) +
        geom_bar(stat = "identity", position = "fill") +
        scale_fill_brewer(palette = "Set3") +
        labs(title = paste("Relative Abundance at", level, "Level"),
             x = "Sample",
             y = "Relative Abundance",
             fill = level) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0("results/figures/", tolower(level), "_composition.pdf"), 
           p, width = 10, height = 6)
}

# Function to plot network
plot_network <- function(network_results) {
    # Extract network and clusters
    net <- network_results$network
    clusters <- network_results$clusters
    
    # Create layout
    layout <- layout_with_fr(net)
    
    # Create plot
    pdf("results/figures/network_plot.pdf", width = 10, height = 10)
    plot(net,
         vertex.color = clusters$membership,
         vertex.size = 5,
         vertex.label = NA,
         edge.width = E(net)$weight,
         layout = layout,
         main = "Microbial Network")
    dev.off()
}

# Function to plot differential abundance
plot_differential_abundance <- function(ps, level) {
    # Convert to DESeq2 format
    deseq_data <- phyloseq_to_deseq2(ps, ~Group)
    
    # Run DESeq2
    deseq_results <- DESeq(deseq_data)
    res <- results(deseq_results)
    
    # Create volcano plot
    volcano_data <- data.frame(
        log2FoldChange = res$log2FoldChange,
        pvalue = res$pvalue,
        padj = res$padj
    )
    
    p <- ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(pvalue))) +
        geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not significant"))) +
        scale_color_manual(values = c("gray", "red")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
        labs(title = paste("Differential Abundance at", level, "Level"),
             x = "Log2 Fold Change",
             y = "-Log10 p-value") +
        theme_minimal()
    
    ggsave(paste0("results/figures/differential_abundance_", tolower(level), ".pdf"),
           p, width = 8, height = 6)
} 