# Load dataset (Assuming it's a CSV file)
Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylation_filtered_means.csv", header = TRUE, row.names = 1)

# Calculate row-wise variance
row_variances <- apply(Methylation, 1, var)

# Filter rows with variance >= 0.3
filtered_dataset <- Methylation[row_variances >= 0.08, ]

# Save the filtered dataset
write.csv(filtered_dataset, "filtered_dataset_Methylation_less than 0.08.csv", row.names = TRUE)

# Print summary
cat("Original rows:", nrow(Methylation), "\n")
cat("Filtered rows:", nrow(filtered_dataset), "\n")

