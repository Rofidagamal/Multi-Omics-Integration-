#  - Multi Omics Course                   # 
#  - Gene Expression Profiling of Prostate#
#  - TCGA-PRAD project - RNAseq data      # 
#  - Assignment 4                         #    
#  - 2024- 1-24                           #
#  - Copyright: @Mohamed Hamed            #
###########################################
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

####load libraries ########  
library(readr)

library("vsn")

library(DESeq2)


##### load the mRNA-Seq data #####
files <- read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\filtered_dataset_RNA_less than 0.3.csv", row.names = 1)

temp <- read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\filtered_dataset_RNA_less than 0.3.csv", row.names = 1)

#########################  Assignment #########################
## load the data using the two other approaches in the following link
# a professional example to follow is here using two ways 
#( reading from files , TCGAbiolink package) >>>  https://www.biostars.org/p/9500223/

### General Reminder
# don't forget to master the tidyvers packages https://www.tidyverse.org/packages/


# check duplciation of of gene symbols?  
x=duplicated(rownames(files))  
sum(x)


###### load the mrna sample sheets  # sample sheets
pheno <- read.csv("C:\\Users\\Electronica Care\\Downloads\\metadata-new (1).csv",header=T)
View(pheno)
library(tidyverse)
metadata <- pheno %>%
  separate(Title, into = c("patient", "sample.type"), sep = ",", fill = "right")


table(metadata$sample.type)

#we will rename the columns of our exp data with the sample ids columns of the pheno file
#however we need to match the file ids to 
file.id= colnames(files)
view((files))
file.ids.metadata=metadata$Accession
index.files=match(file.id,file.ids.metadata)
index.files
colnames(files)=metadata$X [index.files]

save(files,  metadata, file="MODA_RNA_Seq.RDATA")


#### Exploratory analysis + filtration process : plz do it  ####################
# box plot , histogram , PCA
## remove the least 10 % variable genes ### plz do it here


#### Do differential analysis using Deseq2 package as it works on readcounts  ##########
## or use Limma package instead >> you can follow this tutorial here https://www.bioconductor.org/packages/devel/workflows/vignettes/RNAseq123/inst/doc/limmaWorkflow.html
##. read this link : https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

table(metadata$sample.type)

cond1="lesion" 
cond2="NAWM"

ncol(files) == nrow(metadata)
# Remove Zeros
df_clean <- df[apply(df, 1, function(row) all(!is.na(row) & row != 0)), ]

dds = DESeqDataSetFromMatrix( countData = files, colData = metadata , design = ~ sample.type)
dds.run = DESeq(dds)
normalized_counts <- counts(dds.run, normalized=TRUE)

write.csv(normalized_counts,"normalized_RNA.csv")