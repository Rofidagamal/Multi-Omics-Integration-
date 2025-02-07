## ----global_options, include=FALSE----------------------------------------------------------------------------------
library(knitr)
knitr::opts_chunk$set(dpi = 100, echo= TRUE, warning=FALSE, message=FALSE, fig.align = 'center',
                      fig.show=TRUE, fig.keep = 'all', fig.height= 6, fig.width=7)


## -------------------------------------------------------------------------------------------------------------------
library(mixOmics) # import the mixOmics library

set.seed(123) # for reproducibility, remove for normal use

RNA = read.csv('/Users/hageradel/Desktop/Data/Diablo/FinaaalllAnnotation_RNA_mapped_82.csv', header = T, row.names = 1,check.names=F)
Methylation = read.csv('/Users/hageradel/Desktop/Data/Diablo/FinallllAnnotation_Metyl_mapped_8000.csv', header = T, row.names = 1, check.names = F)
metadata= read.csv('/Users/hageradel/Desktop/Data/Diablo/Metadata_RNA_Methyl.csv', header = T)

## -------------------------------------------------------------------------------------------------------------------
# transpose the data (RNA, Methylation)
RNA_t=t(RNA)
Methylation_t=t(Methylation)

# set a list of all the X dataframes
finaldata =list(mRNA=RNA_t,
           Meth= Methylation_t)

# check their dimensions
lapply(data, dim) 

#sepration of title column into 2 columns (Patient", "sample.type)
library(tidyr)
library(dplyr)
metadata_new <- metadata %>%
  separate(Title, into = c("Patient", "sample.type"), sep = ",")

# set the response variable as the Y dataframe
Y = as.factor(metadata_new[1:18,]$sample.type)
summary(Y)

## ---- fig.show = "hold", out.width = "33%", fig.cap = "FIGURE 1: Circle Correlation Plots for pairwise PLS models on the breast TCGA data. Only displays the top 25 features for each dimension, subsetting by those with a correlation above 0.5. "----
# select arbitrary values of features to keep
list.keepX = c(81, 81) # Select 200 features per component
list.keepY = c(3481, 3481) # Keep 200 features from the response variable

str(finaldata[["mRNA"]])   # Check structure
str(finaldata[["Meth"]])

finaldata[["mRNA"]] <- as.matrix(finaldata[["mRNA"]])
finaldata[["Meth"]] <- as.matrix(finaldata[["Meth"]])

sum(is.na(finaldata[["mRNA"]]))  # Count missing values
sum(is.na(finaldata[["Meth"]]))

any(!sapply(finaldata[["mRNA"]], is.numeric))  # Check for non-numeric columns
any(!sapply(finaldata[["Meth"]], is.numeric))

finaldata[["mRNA"]] <- apply(finaldata[["mRNA"]], 2, as.numeric)
finaldata[["Meth"]] <- apply(finaldata[["Meth"]], 2, as.numeric)


pls1 <- spls(finaldata[["mRNA"]], finaldata[["Meth"]], keepX = list.keepX, keepY = list.keepY)

# generate three pairwise PLS models
pls1 <- spls(finaldata[["mRNA"]], finaldata[["Meth"]], keepX = list.keepX, keepY = list.keepY) 

plotVar(pls1, cutoff = 0.5, title = "(a) mRNA vs Meth", legend = c("mRNA", "Meth"), # plot features of first PLS
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

## -------------------------------------------------------------------------------------------------------------------

pls1 <- spls(finaldata[["mRNA"]], finaldata[["Meth"]], ncomp = 18, keepX = 81, keepY = 3481)

cor(pls1$variates$X, pls1$variates$Y) # calculate correlation of miRNA and mRNA
## -------------------------------------------------------------------------------------------------------------------

design = matrix(0.1, ncol = length(finaldata), nrow = length(finaldata), # for square matrix filled with 0.1s
                dimnames = list(names(finaldata), names(finaldata)))
diag(design) = 0 # set diagonal to 0s
design
## -------------------------------------------------------------------------------------------------------------------
basic.diablo.model = block.splsda(X = finaldata, Y = Y, ncomp = 10, design = design) # form basic DIABLO model

## ---- fig.cap = "FIGURE 2: Choosing the number of components in `block.plsda` using `perf()` with 10 × 10-fold CV function in the data. Classification error rates (overall and balanced, see Section 7.3) are represented on the y-axis with respect to the number of components on the x-axis for each prediction distance presented in PLS-DA"----
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', folds = 9, nrepeat = 10) # run component number tuning with repeated CV
plot(perf.diablo) # plot output of tuning

#The most reliable way to determine if the Overall.BER is decreasing is to look at the actual numerical values.
# Or similar, depending on the structure of the perf.diablo object
perf.diablo$error.rate 

## -------------------------------------------------------------------------------------------------------------------
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] # set the optimal ncomp value
perf.diablo$choice.ncomp$WeightedVote # show the optimal choice for ncomp for each dist metric

#save splsda model
save(basic.diablo.model,list.keepX, file = '/Users/hageradel/Desktop/Data/Diablo/splsDA.RData')
## -------------------------------------------------------------------------------------------------------------------
load('/Users/hageradel/Desktop/Data/Diablo/splsDA.RData')
## -------------------------------------------------------------------------------------------------------------------
# ---- eval=FALSE, include = FALSE-----------------------------------------------------------------------------------
test.keepX = list (mRNA = c(5:9, seq(10, 18, 2), seq(20,30,5)),
                 Meth = c(5:9, seq(10, 18, 2), seq(20,30,5)))
     
t1 = proc.time()
tune.plsaDA = tune.block.splsda(X = finaldata, Y = Y, ncomp = 3,
                               test.keepX = test.keepX, design = design,
                               validation = 'Mfold', folds = 9, nrepeat = 1,
                               'cups' = 2, dist = "centroids.dist")
t2 = proc.time()
running_time = t2 - t1; running_time
list.keepX = tune.plsaDA$choice.keepX
list.keepX

save(tune.plsaDA,list.keepX, file = '/Users/hageradel/Desktop/Data/Diablo/splsDA model after test.RData')
load('/Users/hageradel/Desktop/Data/Diablo/splsDA model after test.RData')
## -------------------------------------------------------------------------------------------------------------------

list.keepX = tune.plsaDA$choice.keepX # set the optimal values of features to retain
list.keepX

## -------------------------------------------------------------------------------------------------------------------
final.diablo.model = block.splsda(X = finaldata, Y = Y, ncomp = 3, # set the optimised DIABLO model
                                  keepX = list.keepX, design = design)
final.diablo.model$design # design matrix for the final model

## -------------------------------------------------------------------------------------------------------------------
selectVar(final.diablo.model, block = 'mRNA', comp = 1)$mRNA$name # the features selected to form the first component

## ---- fig.cap = "FIGURE 3: 
#Diagnostic plot from multiblock sPLS-DA applied on the data. 
#Samples are represented based on the specified component (here `ncomp = 1`) for each data set (mRNA, miRNA and protein). 
#Samples are coloured by patenit subtype and 95% confidence ellipse plots are represented."
plotDiablo(final.diablo.model, ncomp = 1)

## ---- fig.cap = "FIGURE 4: 
#Sample plot from multiblock sPLS-DA performed on the data. 
#The samples are plotted according to their scores on the first 2 components for each data set. 
#Samples are coloured by cancer subtype"----
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, title = 'DIABLO Sample Plots')

## ---- fig.cap = "FIGURE 5: 
#Arrow plot from multiblock sPLS-DA performed on the data. 
#The samples are projected into the space spanned by the first two components for each data set then overlaid across data sets.
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, title = 'DIABLO')

## ---- fig.cap = "FIGURE 6: 
#Correlation circle plot from multiblock sPLS-DA performed on the data. 
#Variable types are indicated with different symbols and colours, and are overlaid on the same plot.
plotVar(final.diablo.model, var.names = FALSE, style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2), col = c('darkorchid', 'brown1'))

## ---- fig.cap = "FIGURE 7: 
#Circos plot from multiblock sPLS-DA performed on the data. 
#The plot represents the correlations greater than 0.7 between variables of different types, represented on the side quadrants
circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'lightgreen'),
           color.cor = c("black","yellow"), size.labels = 1.5)

## ---- eval = TRUE, fig.cap = "FIGURE 8: 
#Relevance network for the variables selected by multiblock sPLS-DA performed on the data on component 1. 
#Each node represents a selected with colours indicating their type. The colour of the edges represent positive or negative correlations

png('/Users/hageradel/Desktop/Data/Diablo/Finalll figgg/finalnetwork_plot.png', width = 1000, height = 800) # Adjust width and height in pixels
network(final.diablo.model, blocks = c(1,2), color.node = c('brown1', 'lightgreen'), cutoff = 0.4)
dev.off()

## ---- fig.cap = "FIGURE 9: 
#Loading plot for the variables selected by multiblock sPLS-DA performed on the data on component 1. 
#The most important variables (according to the absolute value of their coefficients) are ordered from bottom to top. 
#As this is a supervised analysis, colours indicate the class for which the median expression value is the highest for each feature
plotLoadings(final.diablo.model, comp = 3, contrib = 'max', method = 'median')

## ---- eval = TRUE, fig.cap = "FIGURE 10: 
#Clustered Image Map for the variables selected by multiblock sPLS-DA performed on the data on component 1. 
#By default, Euclidean distance and Complete linkage methods are used. 
#The CIM represents samples in rows (indicated by their breast cancer subtype on the left hand side of the plot) and selected features in columns (indicated by their data type at the top of the plot).
png('/Users/hageradel/Desktop/Data/Diablo/Finalll figgg/final cimDiablo.png', width = 1000, height = 800) # Adjust width and height in pixels
cimDiablo(final.diablo.model)
dev.off()

## -------------------------------------------------------------------------------------------------------------------
perf.diablo = perf(final.diablo.model, validation = 'Mfold', M = 10, nrepeat = 10, 
                   dist = 'centroids.dist') # run repeated CV performance evaluation

perf.diablo$MajorityVote.error.rate
perf.diablo$WeightedVote.error.rate

## ---- fig.cap = "FIGURE 11: 
#ROC and AUC based on multiblock sPLS-DA performed on the data for the miRNA data set after 2 components. 
#The function calculates the ROC curve and AUC for one class vs. the others."
auc.splsda = auroc(final.diablo.model, roc.block = "miRNA", roc.comp = 2, print = FALSE)

# Install required packages if not installed
if (!require("pROC")) install.packages("pROC", dependencies = TRUE)
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)

# Load Libraries
library(pROC)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Generate Sample Data
n <- 100  # Number of samples
actual <- sample(c(0, 1), n, replace = TRUE)  # Binary labels (0 = Negative, 1 = Positive)
pred1 <- runif(n, min = 0, max = 1)  # Predictions for Model 1
pred2 <- runif(n, min = 0, max = 1)  # Predictions for Model 2

# Compute ROC Curves
roc1 <- roc(actual, pred1)  
roc2 <- roc(actual, pred2)  

# Extract Data for Plotting
roc1_df <- data.frame(TPR = rev(roc1$sensitivities), FPR = rev(1 - roc1$specificities), Model = "mRNA Model")
roc2_df <- data.frame(TPR = rev(roc2$sensitivities), FPR = rev(1 - roc2$specificities), Model = "Meth Model")

# Combine into One Data Frame
roc_data <- rbind(roc1_df, roc2_df)

# Plot ROC Curves
ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1.2) +  # Line for ROC curve
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +  # Diagonal reference line
  annotate("text", x = 0.6, y = 0.2, label = paste("AUC mRNAl 1:", round(auc(roc1), 3)), size = 6, color = "blue") +
  annotate("text", x = 0.6, y = 0.1, label = paste("AUC Meth 2:", round(auc(roc2), 3)), size = 6, color = "green") +
  labs(title = "ROC and AUC in R", x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal() +
  theme(text = element_text(size = 16),
        plot.title = element_text(face = "bold", size = 20),
        legend.position = "bottom")

## -------------------------------------------------------------------------------------------------------------------
predict.diablo = predict(final.diablo.model, newdata = finaldata)

## -------------------------------------------------------------------------------------------------------------------
confusion.mat = get.confusion_matrix(truth = Y,
                                     predicted = predict.diablo$WeightedVote$centroids.dist[,2])
confusion.mat

## -------------------------------------------------------------------------------------------------------------------
get.BER(confusion.mat)

png("heatmap.png", width = 1000, height = 800) # Adjust width and height in pixels
library(ggplot2)
library(reshape2) # For melting the matrix

# Melt the confusion matrix for ggplot2
confusion.mat.melted <- melt(confusion.mat)
names(confusion.mat.melted) <- c("Actual", "Predicted", "Count")

# Create the heatmap
ggplot(data = confusion.mat.melted, aes(x = Predicted, y = Actual, fill = Count)) +
  geom_tile() + # Use tiles for the heatmap
  geom_text(aes(label = Count), color = "black") + # Add count labels
  scale_fill_gradient(low = "white", high = "darkgreen") + # Color gradient
  theme_bw() + # Use a clean theme
  labs(title = "Confusion Matrix", x = "Predicted", y = "Actual") # Add title and labels

dev.off()
## -------------------------------------------------------------------------------------------------------------------
#print the accuracy, specificity, sensitivity and precision of confusion.mat
confusion.mat <- matrix(c(7, 2, 0, 9), nrow = 2, 
                        dimnames = list(Actual = c("lesion", "NAWM"),
                                        Predicted = c("predicted.as.lesion", "predicted.as.NAWM")))

# Calculate metrics
TP_lesion <- confusion.mat["lesion", "predicted.as.lesion"]
TN_lesion <- confusion.mat["NAWM", "predicted.as.NAWM"]  # Corrected TN for lesion
FP_lesion <- confusion.mat["NAWM", "predicted.as.lesion"] # Corrected FP for lesion
FN_lesion <- confusion.mat["lesion", "predicted.as.NAWM"]

TP_NAWM <- confusion.mat["NAWM", "predicted.as.NAWM"]
TN_NAWM <- confusion.mat["lesion", "predicted.as.lesion"] # Corrected TN for NAWM
FP_NAWM <- confusion.mat["lesion", "predicted.as.NAWM"] # Corrected FP for NAWM
FN_NAWM <- confusion.mat["NAWM", "predicted.as.lesion"]

accuracy <- (TP_lesion + TP_NAWM) / sum(confusion.mat)

sensitivity_lesion <- TP_lesion / (TP_lesion + FN_lesion)
specificity_lesion <- TN_lesion / (TN_lesion + FP_lesion)
precision_lesion <- TP_lesion / (TP_lesion + FP_lesion)

sensitivity_NAWM <- TP_NAWM / (TP_NAWM + FN_NAWM)
specificity_NAWM <- TN_NAWM / (TN_NAWM + FP_NAWM)
precision_NAWM <- TP_NAWM / (TP_NAWM + FP_NAWM)


# Print the results
cat("Accuracy:", accuracy, "\n")
cat("Sensitivity (lesion):", sensitivity_lesion, "\n")
cat("Specificity (lesion):", specificity_lesion, "\n")
cat("Precision (lesion):", precision_lesion, "\n")
cat("Sensitivity (NAWM):", sensitivity_NAWM, "\n")
cat("Specificity (NAWM):", specificity_NAWM, "\n")
cat("Precision (NAWM):", precision_NAWM, "\n")

# Or, print in a more organized table format:

results_table <- data.frame(
  Lesion = c(sensitivity_lesion, specificity_lesion, precision_lesion),
  NAWM = c(sensitivity_NAWM, specificity_NAWM, precision_NAWM)
)
rownames(results_table) <- c("Sensitivity", "Specificity", "Precision")

print(results_table)

## -------------------------------------------------------------------------------------------------------------------
selected_vars <- selectVar(final.diablo.model, comp = 3)  # Extract for component 1
selected_vars$mRNA  # Selected features from X dataset
selected_vars$Meth  # Selected features from Y dataset (if applicable)

write.csv(selected_vars$mRNA$name,'/Users/hageradel/Desktop/Data/Diablo/mixOmics Genes.mRNA.csv') 
write.csv(selected_vars$Meth$name,'/Users/hageradel/Desktop/Data/Diablo/mixOmics Genes.Meth.csv')         


