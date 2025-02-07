# Load dataset (Assuming it's a CSV file)
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_norm_counts_FPKM_GRCh38.p13_NCBI.tsv",sep='\t',header = TRUE, row.names = 1,check.names = F)


# Calculate the mean of the entire dataset
overall_mean <- mean(as.matrix(RNA), na.rm = TRUE)

# Filter rows where the row mean is greater than or equal to the overall mean
df_clean <- RNA[rowMeans(RNA, na.rm = TRUE) > 10, ]

write.csv(df_clean,"RNA_filtered_means.csv")



# Load dataset (Assuming it's a CSV file)
Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\filtered_dataset_Methylation_less than 0.08.csv", header = TRUE, row.names = 1)


# Calculate the mean of the entire dataset
overall_mean <- mean(as.matrix(Methylation), na.rm = TRUE)

# Filter rows where the row mean is greater than or equal to the overall mean
df_clean <- Methylation[rowMeans(Methylation, na.rm = TRUE) >0.2, ]

write.csv(df_clean,"Methylation_filtered_means.csv")

