# Load the package
library(org.Hs.eg.db)
# Load your CSV file into R
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Rna_degs.csv", header = TRUE, row.names = 1)
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


write.csv(RNA,"Annotation_RNA_DEGs.csv")

