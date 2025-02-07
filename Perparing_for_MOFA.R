# Load the necessary package
library(org.Hs.eg.db)

# Read the RNA normalized data
RNA_normalized <- read.csv("C:\\Users\\Electronica Care\\Downloads\\normalized_RNAcounts.csv", header = TRUE, row.names = 1)

# Assume Entrez gene IDs are in the row names (already done in the code above)
entrez_ids <- rownames(RNA_normalized)

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID")

# Add gene symbols to the original data
RNA_normalized$GeneSymbol <- gene_symbols

# Remove rows with missing gene symbols (null values)
rna_data_mapped <- RNA_normalized[!is.na(RNA_normalized$GeneSymbol),]

# Calculate the row means
row_means <- rowMeans(rna_data_mapped[, -ncol(rna_data_mapped)], na.rm = TRUE)

# Calculate the number of zeros in each row
num_zeros <- rowSums(rna_data_mapped[, -ncol(rna_data_mapped)] == 0)

# Calculate the proportion of zeros in each row
prop_zeros <- num_zeros / ncol(rna_data_mapped[, -ncol(rna_data_mapped)])

# Set the threshold as 10% of the row mean
threshold <- 0.1 * row_means

# Filter rows where the proportion of zeros is greater than the threshold
filtered_rna_data <- rna_data_mapped[prop_zeros <= threshold, ]

# Output the number of rows before and after filtering
cat("Number of rows before filtering:", nrow(rna_data_mapped), "\n")
cat("Number of rows after filtering:", nrow(filtered_rna_data), "\n")

# Save the filtered data to a new CSV file
write.csv(filtered_rna_data, "Filtered_RNA_standard_normalized.csv", row.names = TRUE)

########################################################################################
library(dplyr)

Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Annotated_Methylation_Data_all_20000.csv", header = TRUE)
Methylation <- Methylation[!is.na(Methylation$UCSC_RefGene_Name),]
Methylation <- Methylation[Methylation$UCSC_RefGene_Name!="",]

library(stringr)
colnames(Methylation) <- str_trim(colnames(Methylation))

colnames(Methylation)[colnames(Methylation) == "UCSC_RefGene_Name"] <- "UCSC_RefGene_Name"
# Aggregate by 'UCSC_RefGene_Name' by calculating the average for duplicate entries
aggregated_data <- Methylation %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(across(.cols = starts_with("GSM"), .fns = mean, na.rm = TRUE))

# Convert to data frame to preserve row names
aggregated_data <- as.data.frame(aggregated_data)

# Set the row names
rownames(aggregated_data) <- aggregated_data$UCSC_RefGene_Name

# Remove the UCSC_RefGene_Name column
aggregated_data <- aggregated_data %>% select(-UCSC_RefGene_Name)

# Calculate the row means
row_means <- rowMeans(aggregated_data[, -ncol(aggregated_data)], na.rm = TRUE)

# Calculate the number of zeros in each row
num_zeros <- rowSums(aggregated_data[, -ncol(aggregated_data)] == 0)

# Calculate the proportion of zeros in each row
prop_zeros <- num_zeros / ncol(aggregated_data[, -ncol(aggregated_data)])

# Set the threshold as 10% of the row mean
threshold <- 0.1 * mean(row_means)

# Filter rows where the proportion of zeros is greater than the threshold
filtered_aggregated_data_methy <- aggregated_data[prop_zeros <= threshold, ]

# Output the number of rows before and after filtering
cat("Number of rows before filtering:", nrow(aggregated_data), "\n")
cat("Number of rows after filtering:", nrow(filtered_aggregated_data_methy), "\n")
filtered_aggregated_data_methy=filtered_aggregated_data_methy[-c(1:10),]

write.csv(filtered_aggregated_data_methy,"Methylated_stan_norm.csv")