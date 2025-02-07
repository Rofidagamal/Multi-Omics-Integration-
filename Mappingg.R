# Load the package
library(org.Hs.eg.db)
# Load your CSV file into R
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\RNA 82 degs.csv", header = TRUE, row.names = 1)
# Assume the Entrez gene IDs are in a column named "EntrezID"
entrez_ids <- rownames(RNA)

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID")

# Add gene symbols to the original data
RNA$GeneSymbol <- gene_symbols


length(!is.na(unique(RNA$GeneSymbol)))

rna_data=RNA$GeneSymbol[!is.na(RNA$GeneSymbol)]
length(rna_data)
length(unique(rna_data))


aggregated_data <- RNA %>%
  group_by(GeneSymbol) %>%
  summarise(across(.cols = starts_with("GSM"), .fns = mean, na.rm = TRUE))


Methyl=read.csv("mapped_data_methyl_8000.csv",header = TRUE, row.names = 1)
aggregated_Methyl <- Methyl %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(across(.cols = starts_with("GSM"), .fns = mean, na.rm = TRUE))

aggregated_Methyl=as.data.frame(aggregated_Methyl)
aggregated_Methyl=aggregated_Methyl[!is.na(aggregated_Methyl$UCSC_RefGene_Name),]
rownames(aggregated_Methyl)=aggregated_Methyl$UCSC_RefGene_Name
aggregated_Methyl <- aggregated_Methyl[, !(names(aggregated_Methyl) %in% "UCSC_RefGene_Name")]
aggregated_Methyl=aggregated_Methyl[-(c(1:8)),]



metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

aggregated_Methy=as.data.frame(aggregated_Methyl)
# Rename columns in dataset if they exist in the metadata
colnames(aggregated_Methy) <- rename_vector[colnames(aggregated_Methy)]
colnames(aggregated_data) <- rename_vector[colnames(aggregated_data)]


aggregated_rna=aggregated_data[,colnames(aggregated_Methy)]

write.csv(aggregated_data,"Annotation_RNA_mapped_82.csv")
write.csv(aggregated_Methy,"Annotation_Metyl_mapped_8000.csv")



###########################################################################################
# Load the package
library(org.Hs.eg.db)
# Load your CSV file into R
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_raw_counts_GRCh38.p13_NCBI (1) (1).tsv", header = TRUE, row.names = 1,sep="\t")
# Assume the Entrez gene IDs are in a column named "EntrezID"
entrez_ids <- rownames(RNA)

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID")

# Add gene symbols to the original data
RNA$GeneSymbol <- gene_symbols


length(!is.na(unique(RNA$GeneSymbol)))

rna_data=RNA$GeneSymbol[!is.na(RNA$GeneSymbol)]
rna_data=RNA[!is.na(RNA$GeneSymbol),]

length(rna_data)
length(unique(rna_data))

library(dplyr)
aggregated_data <- rna_data %>%
  group_by(GeneSymbol) %>%
  summarise(across(.cols = starts_with("GSM"), .fns = mean, na.rm = TRUE))
aggregated_data=as.data.frame(aggregated_data)
rownames(aggregated_data)=aggregated_data$GeneSymbol
aggregated_data <- aggregated_data[, !(names(aggregated_data) %in% "GeneSymbol")]




Methyl=read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\Annotated_Methylation_Data_all_20000.csv",header = TRUE, row.names = 1)
aggregated_Methyl <- Methyl %>%
  group_by(UCSC_RefGene_Name) %>%
  summarise(across(.cols = starts_with("GSM"), .fns = mean, na.rm = TRUE))

aggregated_Methyl=as.data.frame(aggregated_Methyl)
aggregated_Methyl=aggregated_Methyl[!is.na(aggregated_Methyl$UCSC_RefGene_Name),]
rownames(aggregated_Methyl)=aggregated_Methyl$UCSC_RefGene_Name
aggregated_Methyl <- aggregated_Methyl[, !(names(aggregated_Methyl) %in% "UCSC_RefGene_Name")]
aggregated_Methyl=aggregated_Methyl[-(c(1:10)),]

metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

aggregated_Methy=as.data.frame(aggregated_Methyl)
# Rename columns in dataset if they exist in the metadata
colnames(aggregated_Methy) <- rename_vector[colnames(aggregated_Methy)]
colnames(aggregated_data) <- rename_vector[colnames(aggregated_data)]

write.csv(aggregated_data,"All mapped RNA_Data.csv")
write.csv(aggregated_Methy,"All mapped Methyl_20000_Data.csv")


aggregated_rna=aggregated_data[,colnames(aggregated_Methy)]

write.csv(aggregated_data,"Annotation_RNA_mapped_82.csv")
write.csv(aggregated_Methy,"Annotation_Metyl_mapped_8000.csv")






