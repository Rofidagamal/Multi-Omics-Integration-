# ====================================
# Differential Methylation Analysis Script
# Author: [Your Name]
# Description: This script processes methylation data, performs differential analysis, 
# and identifies differentially methylated regions (DMRs) with statistical significance.
# ====================================

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install missing packages if necessary
required_packages <- c("readr", "multtest", "genefilter", "tidyverse", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load libraries
library(readr)
library(multtest)
library(genefilter)
library(tidyverse)
library(ggplot2)

# ====================================
# Function for Multiple Testing Correction
# ====================================
correctPvalueandReturnAll <- function(tt.pval, method = "BH") {
  mt <- mt.rawp2adjp(tt.pval, proc = method)
  adjp <- mt$adjp[order(mt$index), ]
  return(adjp[, 2])  # Adjusted p-values
}

# ====================================
# Load and Preprocess DNA Methylation Data
# ====================================
load_methylation_data <- function(file_path) {
  if (!file.exists(file_path)) stop("Error: Methylation data file not found.")
  meth_data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
  
  # Remove missing or null rows
  meth_data <- meth_data[complete.cases(meth_data), ]
  return(meth_data)
}

# ====================================
# Load and Preprocess Metadata
# ====================================
load_metadata <- function(metadata_path) {
  if (!file.exists(metadata_path)) stop("Error: Metadata file not found.")
  pheno <- read.csv(metadata_path, header = TRUE, stringsAsFactors = FALSE)

  # Split patient and sample type info if needed
  pheno <- pheno %>%
    separate(Title, into = c("Patient", "SampleType"), sep = ",", fill = "right")

  # Clean column names (replace spaces with dots)
  names(pheno) <- gsub(" ", ".", names(pheno))
  
  return(pheno)
}

# ====================================
# Match Sample IDs Between Metadata & Methylation Data
# ====================================
match_samples <- function(methyl_data, metadata) {
  lesion_samples <- metadata$Sample[grepl("lesion", metadata$SampleType)]
  nawm_samples <- metadata$Sample[grepl("NAWM", metadata$SampleType)]
  
  # Extract matched columns
  x <- methyl_data[, lesion_samples, drop = FALSE]
  y <- methyl_data[, nawm_samples, drop = FALSE]
  
  # Combine in lesion vs NAWM order
  meth_combined <- cbind(y, x)
  
  return(meth_combined)
}

# ====================================
# Perform Differential Methylation Analysis
# ====================================
differential_methylation_analysis <- function(methyl_data, case_count, control_count) {
  case_indices <- seq_len(case_count)
  control_indices <- (case_count + 1):(case_count + control_count)
  
  # Compute log fold-change (LFC)
  lfc.diff <- apply(methyl_data, 1, function(x) mean(x[case_indices]) - mean(x[control_indices]))
  
  # Compute t-test p-values
  f <- factor(c(rep(1, length(case_indices)), rep(2, length(control_indices))))
  t.pval <- rowttests(as.matrix(methyl_data), f)$p.value
  
  # Apply multiple testing correction (Benjamini-Hochberg)
  t.pval.adj <- correctPvalueandReturnAll(t.pval, "BH")
  
  # Combine results
  results <- data.frame(LFC = lfc.diff, P_Value = t.pval, Adjusted_P_Value = t.pval.adj, row.names = rownames(methyl_data))
  
  return(results)
}

# ====================================
# Select Differentially Methylated Regions (DMRs)
# ====================================
select_dmrs <- function(results, lfc_threshold = log2(1.5), pval_threshold = 0.05) {
  significant_dmrs <- results[abs(results$LFC) > lfc_threshold & results$P_Value < pval_threshold, ]
  return(significant_dmrs)
}

# ====================================
# Save Results
# ====================================
save_results <- function(dmrs, output_file) {
  write.csv(dmrs, file = output_file, row.names = TRUE)
  cat("Results saved to:", output_file, "\n")
}

# ====================================
# Visualization: Volcano Plot
# ====================================
plot_volcano <- function(results) {
  ggplot(data = results, aes(x = LFC, y = -log10(P_Value))) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Volcano Plot of Differentially Methylated Sites", x = "Log2 Fold Change", y = "-log10 P-Value")
}

# ====================================
# Main Execution
# ====================================
main <- function() {
  # Define file paths
  methylation_file <- "C:\\Users\\user\\Downloads\\multi_omics\\Methylation_20000.csv"
  metadata_file <- "C:\\Users\\user\\Downloads\\Metadata.csv"

  # Load data
  dna_meth <- load_methylation_data(methylation_file)
  pheno_data <- load_metadata(metadata_file)

  # Match samples
  dna_meth_combined <- match_samples(dna_meth, pheno_data)

  # Save processed data
  save(dna_meth_combined, pheno_data, file = "Methy-seq.RDATA")

  # Perform differential methylation analysis
  results <- differential_methylation_analysis(dna_meth_combined, case_count = 9, control_count = 9)

  # Select significant DMRs
  dmrs <- select_dmrs(results)

  # Save results
  save_results(dmrs, "dmrs_genes.csv")

  # Plot results
  plot_volcano(results)
}

# Run the main function
main()
