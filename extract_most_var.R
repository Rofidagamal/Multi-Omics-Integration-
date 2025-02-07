# ===========================================
# Extract Most Variable Genes from Expression Data
# Author: [Your Name]
# Description: This script extracts the top N most variable genes based on standard deviation.
# ===========================================

# Load necessary library
if (!require("dplyr")) install.packages("dplyr")
library(dplyr)

# ===========================================
# Function: Extract Most Variable Genes
# ===========================================
extract_most_variable_genes <- function(expression_data, num_genes = 20000) {
  # Ensure input is a data frame
  if (!is.data.frame(expression_data)) stop("Error: Input must be a data frame.")
  
  # Check if dataset is empty
  if (nrow(expression_data) == 0) stop("Error: Expression data is empty.")
  
  # Compute standard deviation for each gene
  gene_variability <- apply(expression_data, 1, sd, na.rm = TRUE)
  
  # Select the top N most variable genes
  top_genes <- names(sort(gene_variability, decreasing = TRUE))[1:num_genes]
  
  # Subset the data to keep only the selected genes
  filtered_data <- expression_data[top_genes, , drop = FALSE]
  
  # Return the filtered data frame
  return(filtered_data)
}

# ===========================================
# Main Execution
# ===========================================
main <- function() {
  # Define file paths
  input_file <- "C:\\Users\\user\\Downloads\\GSE224455_series_matrix.txt"
  output_file <- "Methylation_20000.csv"
  
  # Check if file exists
  if (!file.exists(input_file)) stop("Error: Input file not found. Please check the file path.")
  
  # Load expression data
  expression_data <- read.csv(input_file, header = TRUE, row.names = 1, check.names = FALSE, sep = "\t")
  
  # Extract most variable genes
  df_filtered <- extract_most_variable_genes(expression_data, num_genes = 20000)
  
  # Save filtered data
  write.csv(df_filtered, output_file, row.names = TRUE)
  cat("Filtered data saved to:", output_file, "\n")
}

# Run the script
main()
