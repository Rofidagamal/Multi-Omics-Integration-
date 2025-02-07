# Load dataset (Assuming it's a CSV file)
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_raw_counts_GRCh38.p13_NCBI (1).tsv", header = TRUE, row.names = 1,sep='\t')

# Calculate row-wise variance
row_variances <- apply(RNA, 1, var)

# Filter rows with variance >= 0.3
filtered_dataset <- RNA[row_variances >= 0.3, ]

# Save the filtered dataset
write.csv(filtered_dataset, "filtered_dataset_RNA_less than 0.3.csv", row.names = TRUE)

# Print summary
cat("Original rows:", nrow(RNA), "\n")
cat("Filtered rows:", nrow(filtered_dataset), "\n")

