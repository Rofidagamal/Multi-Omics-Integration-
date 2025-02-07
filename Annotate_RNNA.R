# Load the package
library(org.Hs.eg.db)
# Load your CSV file into R
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_raw_counts_GRCh38.p13_NCBI (1).tsv", header = TRUE, row.names = 1,sep='\t')
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

uniprot_ids <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "UNIPROT", keytype = "ENTREZID")
RNA$ProteinID <- uniprot_ids

write.csv(RNA,"Annotation_RNA.csv")

rna_data=RNA$ProteinID[!is.na(RNA$ProteinID)]
length(rna_data)
length(unique(rna_data))



# Load the package
library(org.Hs.eg.db)
# Load your CSV file into R
# RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_raw_counts_GRCh38.p13_NCBI (1).tsv", header = TRUE, row.names = 1,sep='\t')
# Assume the Entrez gene IDs are in a column named "EntrezID"
entrez_ids <- rownames(RNA)

# Convert Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID")

# Add gene symbols to the original data
RNA$GeneSymbol <- gene_symbols


# length(!is.na(unique(RNA$GeneSymbol)))

rna_data=RNA$GeneSymbol[!is.na(RNA$GeneSymbol)]
length(rna_data)
length(unique(rna_data))

uniprot_ids <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "UNIPROT", keytype = "ENTREZID")
RNA$ProteinID <- uniprot_ids

write.csv(RNA,"Annotation_RNA_TMM.csv")

rna_data=RNA$ProteinID[!is.na(RNA$ProteinID)]
length(rna_data)
length(unique(rna_data))
