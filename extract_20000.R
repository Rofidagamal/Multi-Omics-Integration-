# ===========================================
# RNA-Seq Exploratory Analysis & PCA
# Author: [Your Name]
# Description: This script loads gene expression data, performs preprocessing,
# exploratory data analysis (EDA), and Principal Component Analysis (PCA).
# ===========================================

# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install missing packages if necessary
required_packages <- c("readr", "ggfortify", "rgl", "scatterplot3d", "matrixStats", 
                       "devtools", "ComplexHeatmap", "circlize", "tidyverse")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

# Load libraries
library(readr)
library(ggfortify)
library(rgl)
library("scatterplot3d")
library(matrixStats)
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)

# ===========================================
# Function: Load Gene Expression Data
# ===========================================
load_expression_data <- function(file_path) {
  if (!file.exists(file_path)) stop("Error: Expression data file not found.")
  
  exp_data <- read_delim(file_path, delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  
  # Extract gene names and convert to matrix
  genes <- unlist(exp_data[,1])
  exp_data <- exp_data[-1]
  exp_data <- as.matrix(exp_data)
  rownames(exp_data) <- genes
  
  return(exp_data)
}

# ===========================================
# Function: Load Metadata
# ===========================================
load_metadata <- function(metadata_path) {
  if (!file.exists(metadata_path)) stop("Error: Metadata file not found.")
  
  metadata <- read_delim(metadata_path, delim = ",", escape_double = FALSE, trim_ws = TRUE)
  
  # Ensure metadata has proper sample labels
  metadata <- metadata %>%
    separate(Title, into = c("Patient", "SampleType"), sep = ",", fill = "right")
  
  return(metadata)
}

# ===========================================
# Function: Exploratory Data Analysis (EDA)
# ===========================================
perform_eda <- function(expression_data) {
  cat("Summary Statistics:\n")
  print(summary(expression_data))
  
  # Boxplot for overall distribution
  boxplot(expression_data, main = "Gene Expression Boxplot", las = 2)
  
  # Density plot for first 18 samples
  dat <- stack(as.data.frame(expression_data[,1:18]))
  ggplot(dat, aes(x = values, fill = ind)) + 
    geom_density(alpha = 0.3) +
    labs(title = "Density Plot of Gene Expression", x = "Expression Level", y = "Density")
}

# ===========================================
# Function: Perform PCA
# ===========================================
perform_pca <- function(expression_data, metadata) {
  # Remove zero-variance genes
  zero_var_genes <- which(apply(expression_data, 1, var) == 0)
  expression_filtered <- expression_data[-zero_var_genes, ]
  
  # PCA computation
  exp_pca <- prcomp(t(expression_filtered), scale = TRUE)
  
  # PCA Plots
  autoplot(exp_pca, data = metadata, colour = 'SampleType', frame = FALSE) +
    labs(title = "PCA of Gene Expression (No Frame)")
  
  autoplot(exp_pca, data = metadata, colour = 'SampleType', frame = TRUE, frame.type = "norm") +
    labs(title = "PCA of Gene Expression (With Normal Frame)")
  
  autoplot(exp_pca, data = metadata, colour = 'Gender', frame = FALSE) +
    labs(title = "PCA Colored by Gender")
  
  return(exp_pca)
}

# ===========================================
# Function: Extract Most Variable Genes
# ===========================================
extract_top_variable_genes <- function(expression_data, num_genes = 20000) {
  # Compute standard deviation
  sds <- rowSds(expression_data)
  
  # Sort genes by variability
  exp_sds_sorted <- data.frame(expression_data, sds)
  exp_sds_sorted <- exp_sds_sorted[order(exp_sds_sorted$sds, decreasing = TRUE), ]
  
  # Select top variable genes
  top_genes <- rownames(exp_sds_sorted)[1:num_genes]
  
  return(expression_data[top_genes, ])
}

# ===========================================
# Main Execution
# ===========================================
main <- function() {
  # Define file paths
  expression_file <- "C:/Users/user/Downloads/GSE224377_raw_counts_GRCh38.p13_NCBI.tsv"
  metadata_file <- "C:/Users/user/Downloads/Metadata.csv"
  
  # Load data
  exp_data <- load_expression_data(expression_file)
  metadata <- load_metadata(metadata_file)
  
  # Save raw data
  save(exp_data, metadata, file = "MODA_day1.RDATA")
  
  # Perform EDA
  perform_eda(exp_data)
  
  # Perform PCA
  exp_pca <- perform_pca(exp_data, metadata)
  
  # Extract top variable genes
  exp_top20000 <- extract_top_variable_genes(exp_data, num_genes = 20000)
  
  # Save filtered data
  save(exp_top20000, file = "Top20000Genes.RDATA")
  write.csv(exp_top20000, "Top20000Genes.csv", row.names = TRUE)
  
  cat("Data processing completed successfully. Filtered expression data saved.\n")
}

# Run the main function
main()
