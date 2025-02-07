# ============================
# CpG Site Mapping to Gene Symbols
# Author: [Your Name]
# Description: This script maps CpG sites to gene symbols using Bioconductor.
# Dependencies: Bioconductor, minfi, IlluminaHumanMethylationEPICanno.ilmn12.hg19
# ============================

# Install and Load Required Libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install required Bioconductor packages if not installed
if (!require("minfi")) BiocManager::install("minfi")
if (!require("IlluminaHumanMethylationEPICanno.ilmn12.hg19")) {
    BiocManager::install("IlluminaHumanMethylationEPICanno.ilmn12.hg19")
}

# Load necessary libraries
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilmn12.hg19)

# ============================
# Load CpG Site Data
# ============================

# Set file path (Modify as needed)
input_file <- "C:\\Users\\user\\Downloads\\multi_omics\\dmrs_genes.csv"

# Check if the file exists before loading
if (!file.exists(input_file)) {
    stop("Error: Input file not found. Please check the file path.")
}

# Read CpG site data (Ensure the first column contains CpG IDs)
cpg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# Check the structure of the data
if (!"X" %in% colnames(cpg_data)) {
    stop("Error: The input file must have a column named 'X' containing CpG site IDs.")
}

# ============================
# Annotate CpG Sites
# ============================

# Extract CpG site IDs from the dataset
cpg_sites <- cpg_data$X

# Load annotation data for Illumina EPIC array
annotation_data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilmn12.hg19)

# Perform mapping by matching CpG sites
mapped_cpg <- annotation_data[rownames(annotation_data) %in% cpg_sites, ]

# Check if any mappings were found
if (nrow(mapped_cpg) == 0) {
    stop("Error: No matching CpG sites found in the annotation dataset.")
}

# Convert to a data frame for easy export
mapped_cpg_df <- data.frame(
    CpG_ID = rownames(mapped_cpg),
    Gene_Symbol = mapped_cpg$UCSC_RefGene_Name,  # UCSC RefGene column contains gene symbols
    Gene_Type = mapped_cpg$UCSC_RefGene_Group  # Gene region (e.g., TSS200, Body, etc.)
)

# ============================
# Save Mapped Data
# ============================

# Define output file path
output_file <- "mapped_cpg_sites.csv"

# Save mapped CpG sites to a CSV file
write.csv(mapped_cpg_df, output_file, row.names = FALSE)

# Success message
cat("CpG site mapping completed successfully. Output saved to:", output_file, "\n")
