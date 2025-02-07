library(MatrixGenerics)
library(matrixStats)
library(IRanges)
library(GenomeInfoDb)

library(dynamicTreeCut)
library(fastcluster)
library(WGCNA)
# Install the corrplot package if it's not already installed
install.packages("corrplot")

# Load the corrplot package
library(corrplot)
install.packages("remotes")
remotes::install_github("kevinblighe/CorLevelPlot")

install.packages("WGCNA")
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(ggplot2)
library(dplyr)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(RColorBrewer)

# 1.RNA-Seq data
RNA <- read.csv("C:\\Users\\user\\Downloads\\multi omics\\GSE224377_raw_counts_GRCh38.p13_NCBI.tsv\\GSE224377_raw_counts_GRCh38.p13_NCBI.tsv", row.names=1, header=T, sep="\t") # Rows are genes, columns are samples

# get metadata 
 
geo_id <- "GSE224377"
gse <- getGEO(geo_id, GSEMatrix= TRUE)
metadata <- pData(phenoData(gse[[1]]))
metadata <- metadata[,c(1,2,21,41,45)]
phenoData <- metadata %>%
  separate(title, into = c("patient", "sample.type"), sep = ",", fill = "right")

# 2. QC-outlier detection 
# detect outlier genes

RNA_t <- goodSamplesGenes(t(RNA))
summary(RNA_t)
RNA_t$allOK  #FALSE means samples or genes have outliers 
table(RNA_t$goodGenes)

#remove genes that are detected as ouliers 
RNA <- RNA[RNA_t$goodGenes==TRUE,]

#detect outlier samples - HCA - method 1
htree <- hclust(dist(t(RNA)), method= "average")
plot(htree)

# pca - method 2
pca <- prcomp(t(RNA)) 
pca.data <- pca$x #information about principal confidence is stored in a slot called 'x'

pca.var <- pca$sdev^2 #variance explained by each principal confidence
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits=2)

pca.data <- as.data.frame(pca.data)
dev.off()
library(PerformanceAnalytics)
ggplot(pca.data, aes(PC1,PC2))+
  geom_point()+
  geom_text(label=rownames(pca.data))+
  labs(x= paste0('PC1: ', pca.var.percent[1],'%'),
       y= paste0('PC2: ', pca.var.percent[2],'%'))

###NOTE: If there are batch effects observed, correct for them before moving ahead

#exclude outlier samples 
samples.to.be.excluded <- c('GSM7021450','GSM7021456')
data.subset <- RNA[,!(colnames(RNA) %in% samples.to.be.excluded)]

# 3. Normalization.....................
#create a deseq2 dataset 
#exclude outlier samples 
library(dplyr)

colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


#fixing column names in colData

names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s','_', names(colData))

#making the rownames and columns identical 
all(rownames(colData) %in% colnames(data.subset)) #same_samples.names
all(rownames(colData) == colnames(data.subset)) #same_sample.order

#create dds
dds <- DESeqDataSetFromMatrix(countData= data.subset,
                              colData= colData,
                              design= ~ sample.type)

## remove all genes with counts <15 in more than 75% of samples (16*0.75=12)
## 16 is the number of variables in our data set now   
## suggested by WGCNA on RNASeq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 12,]
nrow(dds75) # 11236 genes

# perform variance stablization

dds_norm <- vst(dds75)

# get normalized counts 
norm.counts <- assay(dds_norm) %>%
  t()  #t for transpose to be used for WGCNA analysis functions 

# 4. Network Construction 
#choose a set of soft-thresholding powers 
power <- c(c(1:10), seq(from=12, to=50, by=2))

# Call the network topology analysis function

sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)

sft.data <- sft$fitIndices #we are concerned with the SFT.R.sq (maximum r square) and mean.k. (minimum mean)

# visualization to pick power


a1 <- ggplot(sft.data, aes(Power,SFT.R.sq, label= Power))+
  geom_point()+
  geom_text(nudge_y= 0.1)+
  geom_hline(yintercept = 0.8, color= 'red')+
  labs(x= 'Power', y= 'Scale free topology model fit, signed R^2')+
  theme_classic()
  
a2 <- ggplot(sft.data, aes(Power,mean.k., label= Power))+
  geom_point()+
  geom_text(nudge_y= 0.1)+
  geom_hline(yintercept = 0.8, color= 'blue')+
  labs(x= 'Power', y= 'Mean Connectivity')+
  theme_classic()

grid.arrange(a1,a2, nrow=2)

# convert matrix to numeric 
norm.counts[] <- sapply(norm.counts, as.numeric)
soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize 
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power= soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed= 1234,
                 verbose= 3)

cor <- temp_cor

# 5. Module Eigengenes ---------------------
module_eigengenes <- bwnet$MEs

# Print out a preview 
head(module_eigengenes) #samples as the row names and for each sample the module eigengenes have been calculated 

# get number of genes for each module 
table(bwnet$colors)

# Plot the dendrogram and the module colors before and after merging underneath 

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors,bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

# grey module = all genes that doesn't fail into other modules were assigned to the grey module






# 6A. Relate modules to traits ---------------------------------
# module trait associations 




# create traits file - binarize categorical variables 
traits <- colData %>% 
  mutate(sample.type_bin= ifelse(grepl('lesion', sample.type),1,0))%>%
  select(7)

# binarize categorical variables 

colData$seks <- factor(colData$seks, levels = c("F","M"))

seks_out <- binarizeCategoricalColumns(colData$seks,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount= 1)

traits <- cbind(traits,seks_out)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr,nSamples) 
# to identify the modules that are significantly associated with the disease or the gender

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>% 
  column_to_rownames(var='Row.names')

CorLevelPlot(heatmap.data,
             x= names(heatmap.data)[18:19],
             y= names(heatmap.data)[1:17],
             col=c("skyblue","pink"))

# extract the GeneIDs that are in the significant modules 
module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping$Gene_ID = rownames(module.gene.mapping)

red_df = module.gene.mapping[module.gene.mapping$`bwnet$colors` == "red", "Gene_ID"]
yellow_df = module.gene.mapping[module.gene.mapping$`bwnet$colors` == "yellow", "Gene_ID"]

#significant modules associated sample.type by correlation method

# 6B. Intramodular analysis: Identifying driver genes ----------------------
# calculate the module membership and the associated p-values 

# The module membership/intramodular connectivity (hub) is calculated as the correlation of the eigengene 
# This quantifies the similarity if all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

module.membership.measure.pvals[1:10,1:10]

# Calculate the gene significance and associated p-values 

gene.signf.corr <- cor(norm.counts, traits$sample.type_bin, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr.pvals %>% 
  as.data.frame()%>%
  arrange(V1) %>%
  head(25)

# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

rownames(dds)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi") # Required for some functions

library(org.Hs.eg.db)
library(AnnotationDbi)



