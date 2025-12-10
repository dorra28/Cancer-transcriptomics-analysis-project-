#!/usr/bin/env Rscript
#' DESeq2 Differential Expression Analysis for Cancer Transcriptomics
#' 
#' Implements robust differential expression analysis using DESeq2
#' with optimized parameters for cancer RNA-seq data
#' 
#' Reference:
#'   Love, M.I., Huber, W., Anders, S. (2014)
#'   Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2
#'   Genome Biology, 15(12), 550
#' 
#' @author Dorra Rjaibi
#' @date 2025

# Required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(BiocParallel)
})


#' Run DESeq2 Differential Expression Analysis
#' 
#' @param counts_file Path to count matrix file (genes x samples)
#' @param metadata_file Path to clinical metadata file
#' @param design_formula Design formula for DESeq2 (e.g., "~ condition")
#' @param contrast Vector specifying contrast (e.g., c("condition", "tumor", "normal"))
#' @param fdr_threshold FDR cutoff for significance (default: 0.05)
#' @param log2fc_threshold Log2 fold change threshold (default: 1.5)
#' @param output_dir Directory to save results
#' @param n_cores Number of cores for parallel processing
#' 
#' @return DESeq2 results object
run_deseq2_analysis <- function(counts_file,
                                metadata_file,
                                design_formula = "~ condition",
                                contrast = c("condition", "tumor", "normal"),
                                fdr_threshold = 0.05,
                                log2fc_threshold = 1.5,
                                output_dir = "results/deseq2",
                                n_cores = 4) {
  
  cat("\n=== DESeq2 Differential Expression Analysis ===\n")
  cat(sprintf("Date: %s\n", Sys.time()))
  cat(sprintf("DESeq2 version: %s\n", packageVersion("DESeq2")))
  
  # Create output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Setup parallel processing
  register(MulticoreParam(n_cores))
  
  # ============================================================================
  # STEP 1: Load Data
  # ============================================================================
  cat("\n[1] Loading data...\n")
  
  # Load count matrix
  counts_data <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  cat(sprintf("  Loaded %d genes x %d samples\n", nrow(counts_data), ncol(counts_data)))
  
  # Load metadata
  metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
  cat(sprintf("  Loaded metadata for %d samples\n", nrow(metadata)))
  
  # Ensure sample order matches
  metadata <- metadata[colnames(counts_data), , drop = FALSE]
  
  # Verify condition column exists
  condition_col <- contrast[1]
  if (!condition_col %in% colnames(metadata)) {
    stop(sprintf("Condition column '%s' not found in metadata", condition_col))
  }
  
  # Convert to factor with proper reference level
  metadata[[condition_col]] <- factor(metadata[[condition_col]])
  if (contrast[3] %in% levels(metadata[[condition_col]])) {
    metadata[[condition_col]] <- relevel(metadata[[condition_col]], ref = contrast[3])
  }
  
  cat(sprintf("  Condition groups: %s\n", 
              paste(levels(metadata[[condition_col]]), collapse = ", ")))
  
  # ============================================================================
  # STEP 2: Create DESeq2 Dataset
  # ============================================================================
  cat("\n[2] Creating DESeq2 dataset...\n")
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts_data),
    colData = metadata,
    design = as.formula(design_formula)
  )
  
  # Pre-filtering: remove genes with very low counts
  # Keep genes with at least 10 reads in at least 3 samples
  keep <- rowSums(counts(dds) >= 10) >= 3
  dds <- dds[keep, ]
  cat(sprintf("  Retained %d genes after low-count filtering\n", nrow(dds)))
  
  # ============================================================================
  # STEP 3: Run DESeq2 Analysis
  # ============================================================================
  cat("\n[3] Running DESeq2 differential expression analysis...\n")
  cat("  This may take several minutes...\n")
  
  # Run DESeq2 pipeline
  # Uses Wald test for hypothesis testing
  # Applies independent filtering to increase power
  dds <- DESeq(dds, parallel = TRUE)
  
  # Extract results
  res <- results(
    dds,
    contrast = contrast,
    alpha = fdr_threshold,
    parallel = TRUE
  )
  
  # Order by adjusted p-value
  res <- res[order(res$padj), ]
  
  # Summary statistics
  cat("\n  DESeq2 Analysis Summary:\n")
  cat(sprintf("    Total genes tested: %d\n", sum(!is.na(res$padj))))
  cat(sprintf("    Significant genes (FDR < %.2f): %d\n", 
              fdr_threshold, sum(res$padj < fdr_threshold, na.rm = TRUE)))
  cat(sprintf("    Upregulated (log2FC > %.2f): %d\n",
              log2fc_threshold,
              sum(res$log2FoldChange > log2fc_threshold & 
                  res$padj < fdr_threshold, na.rm = TRUE)))
  cat(sprintf("    Downregulated (log2FC < -%.2f): %d\n",
              log2fc_threshold,
              sum(res$log2FoldChange < -log2fc_threshold & 
                  res$padj < fdr_threshold, na.rm = TRUE)))
  
  # ============================================================================
  # STEP 4: Save Results
  # ============================================================================
  cat("\n[4] Saving results...\n")
  
  # Convert results to data frame
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Add significance flags
  res_df$significant <- res_df$padj < fdr_threshold
  res_df$direction <- ifelse(
    res_df$significant,
    ifelse(res_df$log2FoldChange > log2fc_threshold, "Up",
           ifelse(res_df$log2FoldChange < -log2fc_threshold, "Down", "None")),
    "None"
  )
  
  # Save full results
  write.csv(res_df, 
            file.path(output_dir, "deseq2_results_full.csv"),
            row.names = FALSE)
  
  # Save significant DEGs only
  sig_genes <- res_df[res_df$significant & 
                      abs(res_df$log2FoldChange) > log2fc_threshold, ]
  write.csv(sig_genes,
            file.path(output_dir, "deseq2_significant_genes.csv"),
            row.names = FALSE)
  
  cat(sprintf("  Saved results to: %s\n", output_dir))
  
  # ============================================================================
  # STEP 5: Quality Control Plots
  # ============================================================================
  cat("\n[5] Generating quality control plots...\n")
  
  # MA plot (log ratio vs mean)
  png(file.path(output_dir, "ma_plot.png"), width = 800, height = 600, res = 100)
  plotMA(res, ylim = c(-5, 5), main = "MA Plot: DESeq2 Results")
  abline(h = c(-log2fc_threshold, log2fc_threshold), col = "red", lty = 2)
  dev.off()
  
  # Dispersion plot
  png(file.path(output_dir, "dispersion_plot.png"), width = 800, height = 600, res = 100)
  plotDispEsts(dds, main = "Dispersion Estimates")
  dev.off()
  
  # P-value histogram
  png(file.path(output_dir, "pvalue_histogram.png"), width = 800, height = 600, res = 100)
  hist(res$pvalue[res$baseMean > 1], 
       breaks = 50, 
       col = "skyblue",
       main = "P-value Distribution",
       xlab = "P-value",
       ylab = "Frequency")
  dev.off()
  
  # Volcano plot
  generate_volcano_plot(res_df, output_dir, fdr_threshold, log2fc_threshold)
  
  # ============================================================================
  # STEP 6: Normalized Counts and Heatmap
  # ============================================================================
  cat("\n[6] Generating expression heatmaps...\n")
  
  # Get variance-stabilized transformed data
  vsd <- vst(dds, blind = FALSE)
  
  # Save normalized counts
  normalized_counts <- assay(vsd)
  write.csv(normalized_counts,
            file.path(output_dir, "vst_normalized_counts.csv"))
  
  # Heatmap of top DEGs
  if (nrow(sig_genes) > 0) {
    top_genes <- head(sig_genes$gene, min(50, nrow(sig_genes)))
    generate_deg_heatmap(vsd, top_genes, metadata, output_dir)
  }
  
  # PCA plot
  generate_pca_plot(vsd, metadata, condition_col, output_dir)
  
  # ============================================================================
  # STEP 7: Save DESeq2 Object
  # ============================================================================
  cat("\n[7] Saving DESeq2 object...\n")
  
  saveRDS(dds, file.path(output_dir, "deseq2_dds_object.rds"))
  saveRDS(res, file.path(output_dir, "deseq2_results_object.rds"))
  
  cat("\n=== DESeq2 Analysis Complete ===\n")
  cat(sprintf("Results saved to: %s\n", output_dir))
  
  return(list(dds = dds, results = res, results_df = res_df))
}


#' Generate Volcano Plot
#' 
#' @param res_df DESeq2 results data frame
#' @param output_dir Output directory
#' @param fdr_threshold FDR cutoff
#' @param log2fc_threshold Log2 fold change threshold
generate_volcano_plot <- function(res_df, output_dir, fdr_threshold, log2fc_threshold) {
  
  # Prepare data for plotting
  plot_data <- res_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      neg_log10_padj = -log10(padj),
      color_group = case_when(
        padj < fdr_threshold & log2FoldChange > log2fc_threshold ~ "Upregulated",
        padj < fdr_threshold & log2FoldChange < -log2fc_threshold ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  # Create volcano plot
  p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj, color = color_group)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(
      values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey"),
      name = "Differential Expression"
    ) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), 
               linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(fdr_threshold), 
               linetype = "dashed", color = "black", alpha = 0.5) +
    labs(
      title = "Volcano Plot: Differential Gene Expression",
      subtitle = sprintf("FDR < %.2f, |log2FC| > %.2f", fdr_threshold, log2fc_threshold),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  # Add gene labels for top genes
  top_up <- plot_data %>%
    filter(color_group == "Upregulated") %>%
    arrange(desc(log2FoldChange)) %>%
    head(10)
  
  top_down <- plot_data %>%
    filter(color_group == "Downregulated") %>%
    arrange(log2FoldChange) %>%
    head(10)
  
  top_genes <- rbind(top_up, top_down)
  
  if (nrow(top_genes) > 0) {
    p <- p + ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = gene),
      size = 3,
      max.overlaps = 10
    )
  }
  
  ggsave(file.path(output_dir, "volcano_plot.png"), 
         plot = p, width = 10, height = 8, dpi = 300)
}


#' Generate DEG Heatmap
#' 
#' @param vsd Variance-stabilized data
#' @param genes Genes to plot
#' @param metadata Sample metadata
#' @param output_dir Output directory
generate_deg_heatmap <- function(vsd, genes, metadata, output_dir) {
  
  # Extract expression data for selected genes
  mat <- assay(vsd)[genes, ]
  
  # Scale by row (Z-score)
  mat_scaled <- t(scale(t(mat)))
  
  # Prepare annotation
  annotation_col <- data.frame(
    Condition = metadata[[names(metadata)[1]]],
    row.names = colnames(mat)
  )
  
  # Color scheme
  colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
  
  # Generate heatmap
  png(file.path(output_dir, "deg_heatmap.png"), 
      width = 1200, height = 1000, res = 150)
  
  pheatmap(
    mat_scaled,
    color = colors,
    scale = "none",
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    annotation_col = annotation_col,
    show_rownames = ifelse(length(genes) <= 50, TRUE, FALSE),
    show_colnames = FALSE,
    main = "Heatmap: Top Differentially Expressed Genes",
    fontsize = 10,
    fontsize_row = 8
  )
  
  dev.off()
}


#' Generate PCA Plot
#' 
#' @param vsd Variance-stabilized data
#' @param metadata Sample metadata
#' @param condition_col Name of condition column
#' @param output_dir Output directory
generate_pca_plot <- function(vsd, metadata, condition_col, output_dir) {
  
  # Perform PCA
  pca_data <- plotPCA(vsd, intgroup = condition_col, returnData = TRUE)
  percent_var <- round(100 * attr(pca_data, "percentVar"))
  
  # Create PCA plot
  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = get(condition_col))) +
    geom_point(size = 3, alpha = 0.8) +
    labs(
      title = "PCA Plot: Sample Clustering",
      x = sprintf("PC1: %d%% variance", percent_var[1]),
      y = sprintf("PC2: %d%% variance", percent_var[2]),
      color = "Condition"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, "pca_plot.png"), 
         plot = p, width = 8, height = 6, dpi = 300)
}


#' Main execution function
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2) {
    cat("\nUsage: Rscript deseq2_analysis.R <counts_file> <metadata_file> [options]\n")
    cat("\nRequired arguments:\n")
    cat("  counts_file    : Path to count matrix CSV file\n")
    cat("  metadata_file  : Path to clinical metadata CSV file\n")
    cat("\nOptional arguments:\n")
    cat("  --output       : Output directory (default: results/deseq2)\n")
    cat("  --fdr          : FDR threshold (default: 0.05)\n")
    cat("  --log2fc       : Log2 fold change threshold (default: 1.5)\n")
    cat("  --cores        : Number of CPU cores (default: 4)\n")
    cat("\nExample:\n")
    cat("  Rscript deseq2_analysis.R counts.csv metadata.csv --output results/brca --fdr 0.01\n\n")
    quit(status = 1)
  }
  
  # Default parameters
  counts_file <- args[1]
  metadata_file <- args[2]
  output_dir <- "results/deseq2"
  fdr_threshold <- 0.05
  log2fc_threshold <- 1.5
  n_cores <- 4
  
  # Parse optional arguments
  if (length(args) > 2) {
    for (i in 3:length(args)) {
      if (args[i] == "--output" && i < length(args)) {
        output_dir <- args[i + 1]
      } else if (args[i] == "--fdr" && i < length(args)) {
        fdr_threshold <- as.numeric(args[i + 1])
      } else if (args[i] == "--log2fc" && i < length(args)) {
        log2fc_threshold <- as.numeric(args[i + 1])
      } else if (args[i] == "--cores" && i < length(args)) {
        n_cores <- as.integer(args[i + 1])
      }
    }
  }
  
  # Run analysis
  results <- run_deseq2_analysis(
    counts_file = counts_file,
    metadata_file = metadata_file,
    design_formula = "~ condition",
    contrast = c("condition", "tumor", "normal"),
    fdr_threshold = fdr_threshold,
    log2fc_threshold = log2fc_threshold,
    output_dir = output_dir,
    n_cores = n_cores
  )
}


# Execute if run as script
if (!interactive()) {
  main()
}
