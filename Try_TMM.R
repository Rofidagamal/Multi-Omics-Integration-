library("edgeR")
RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_raw_counts_GRCh38.p13_NCBI (1).tsv", header = TRUE, row.names = 1,sep='\t')#mRNA_quantile_normalization
mRNA_normalization <- function(mRNA) {
  #' @param mRNA data_frame, include_mRNA_in_a_data_frame
  dataset_finalM <- as.matrix(mRNA)
  datasetA <- DGEList(counts=dataset_finalM)
  keep.exprs <- filterByExpr(datasetA)
  datasetAfiltered <- datasetA[keep.exprs,keep.lib.sizes=FALSE]
  datasetAnorm <- calcNormFactors(datasetAfiltered, method = "TMM")  
  mRNA_norm <- cpm(datasetAnorm, log=TRUE)
  boxplot(mRNA_norm, las=2, main="")
  title(main="Normalised data",ylab="Log-cpm")
  return(mRNA_norm)
}

#obtain_the_data_into_matrix
RNA = as.matrix(RNA)
#filter based on gene variance
variance_list_of_genes = c()
#calculate_variance_of_genes
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

RNA=mRNA_normalization(RNA)
variance_of_genes= RowVar(RNA)

#save_eliminated_genes_index
eliminated_genes_indices = c()
#Loop_to_eleminate_gene_with_0.3 or < 0.3 low_var_genes
for(i in 1:length(variance_of_genes))
{
  if(variance_of_genes[i]<=0.3)
  {
    eliminated_genes_indices =c(eliminated_genes_indices,i)
  }
  
}
#remove_low_var_genes
filtered_based_on_variance = RNA[-eliminated_genes_indices,]
dim(filtered_based_on_variance)
#filter based on total gene count per row
keep = rowSums((filtered_based_on_variance))>=10
#select_high_var_genes_only
gene_count_criteria_keep = filtered_based_on_variance[keep,]
#select_gene_count_col_only
colnames(gene_count_criteria_keep)=substring(colnames(gene_count_criteria_keep),first=1,last=34)
RNA= gene_count_criteria_keep
#convert_each_count_value_to_integer_number
genes=row.names(RNA)
RNA=apply(RNA,2,as.integer)
row.names(RNA)=genes
RNA=as.data.frame(RNA)