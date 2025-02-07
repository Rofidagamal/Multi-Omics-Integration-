# install.packages("devtools")
# devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

library(EPIC)
rna_data <- read.csv("All mapped RNA_Data.csv", row.names =  1,check.names = F)
result <- EPIC(rna_data, reference="TRef")

cell_fractions <- result$cellFractions
head(cell_fractions)
write.csv(cell_fractions, "EPIC_cell_fractions.csv")

library(ggplot2)
library(reshape2)

cell_fractions=as.data.frame(cell_fractions)
cell_fractions$Sample=rownames(cell_fractions)


# Load required libraries
library(ggplot2)
library(reshape2)

# Reshape data: Convert wide format (samples as rows) to long format (samples as a column)
df <- melt(as.data.frame(cell_fractions), id.vars = "Sample", variable.name = "CellType", value.name = "Proportion")

# Check the reshaped data
head(df)

# Plot cell type proportions per sample
ggplot(df, aes(x=Sample, y=Proportion, fill=CellType)) + 
  geom_bar(stat="identity", position="stack") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title="Cell Type Proportions", x="Sample", y="Proportion", fill="Cell Type")

#######################################################################################
# install.packages("devtools")
# devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)

library(EPIC)
methylation_data <- read.csv("All mapped Methyl_20000_Data.csv", row.names = 1,check.names = F)
result <- EPIC(methylation_data, reference="TRef")

cell_fractions <- result$cellFractions
head(cell_fractions)
write.csv(cell_fractions, "EPIC_cell_fractions_methyl.csv")

library(ggplot2)
library(reshape2)

cell_fractions=as.data.frame(cell_fractions)
cell_fractions$Sample=rownames(cell_fractions)


# Load required libraries
library(ggplot2)
library(reshape2)

# Reshape data: Convert wide format (samples as rows) to long format (samples as a column)
df <- melt(as.data.frame(cell_fractions), id.vars = "Sample", variable.name = "CellType", value.name = "Proportion")

# Check the reshaped data
head(df)

# Plot cell type proportions per sample
ggplot(df, aes(x=Sample, y=Proportion, fill=CellType)) + 
  geom_bar(stat="identity", position="stack") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  labs(title="Cell Type Proportions", x="Sample", y="Proportion", fill="Cell Type")
