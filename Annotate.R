Methylation_data=read.csv("C:\\Users\\Electronica Care\\Downloads\\Annotated_Methylation_Data_DEMs (1).csv", header = TRUE, row.names = 1)
# Install required packages for methylation data analysis
# BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "missMethyl"))

# Load required libraries
# library(minfi)  # For methylation data processing
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # Annotation package for EPIC array
library(missMethyl)  # For downstream analysis like gene set enrichment

# Set CpG sites as row names
# rownames(Methylation_data) <- Methylation_data[[1]]  # Assign first column as row names

# Convert to matrix
methylation_matrix <- as.matrix(Methylation_data)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Select relevant columns (CpG ID & Gene Name)
anno_subset <- anno[, c("Name", "UCSC_RefGene_Name")]

# Convert to a dataframe
anno_df <- as.data.frame(anno_subset)

# Merge methylation data with CpG-to-Gene annotation
mapped_data <- merge(anno_df, as.data.frame(methylation_matrix), 
                     by.x = "Name", by.y = "row.names", all.y = TRUE)

# Save the annotated file
write.csv(mapped_data, "Annotated_Methylation_Data_all_20000.csv")

#######################################################################################
# BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "missMethyl"))

# Install necessary Bioconductor packages if not already installed

Methylation_data=read.csv("C:\\Users\\Electronica Care\\Downloads\\dmrs_genes_count_data.csv", header = TRUE, row.names = 1)
# Install required packages for methylation data analysis
# BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "missMethyl"))

# Load required libraries
library(minfi)  # For methylation data processing
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # Annotation package for EPIC array
library(missMethyl)  # For downstream analysis like gene set enrichment

# Set CpG sites as row names
# rownames(Methylation_data) <- Methylation_data[[1]]  # Assign first column as row names

# Convert to matrix
methylation_matrix <- as.matrix(Methylation_data)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Select relevant columns (CpG ID & Gene Name)
anno_subset <- anno[, c("Name", "UCSC_RefGene_Name")]

# Convert to a dataframe
anno_df <- as.data.frame(anno_subset)

# Merge methylation data with CpG-to-Gene annotation
mapped_data <- merge(anno_df, as.data.frame(methylation_matrix), 
                     by.x = "Name", by.y = "row.names", all.y = TRUE)

# Save the annotated file
write.csv(mapped_data, "Annotated_Methylation_Data_DEMs.csv", row.names = FALSE)


#######################################################################################
# BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "missMethyl"))

# Install necessary Bioconductor packages if not already installed

Methylation_data=read.csv("C:\\Users\\Electronica Care\\Downloads\\dmrs_genes_count_data (1).csv", header = TRUE, row.names = 1)
# Install required packages for methylation data analysis
# BiocManager::install(c("minfi", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "missMethyl"))

# Load required libraries
library(minfi)  # For methylation data processing
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)  # Annotation package for EPIC array
library(missMethyl)  # For downstream analysis like gene set enrichment

# Set CpG sites as row names
# rownames(Methylation_data) <- Methylation_data[[1]]  # Assign first column as row names

# Convert to matrix
methylation_matrix <- as.matrix(Methylation_data)

# Load annotation data
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Select relevant columns (CpG ID & Gene Name)
anno_subset <- anno[, c("Name", "UCSC_RefGene_Name")]

# Convert to a dataframe
anno_df <- as.data.frame(anno_subset)

# Merge methylation data with CpG-to-Gene annotation
mapped_data <- merge(anno_df, as.data.frame(methylation_matrix), 
                     by.x = "Name", by.y = "row.names", all.y = TRUE)

# Save the annotated file
write.csv(mapped_data, "C:/Users/Electronica Care/Desktop/Annotated_Methylation_Data_DEMs.csv")