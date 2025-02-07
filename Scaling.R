# RNA=as.matrix(RNA)
# # 1. Standard Z-Scaling (Z-score normalization)
# mu <- mean(RNA)
# sigma <- sd(RNA)
# z_scaled_data <- (RNA - mu) / sigma
# 





# # 2. Z-Scaling and then Min-Max Scaling to -1 to 1
# 
# # (Z-scaling is already done above, so we just reuse z_scaled_data)
# 
# min_z <- min(z_scaled_data)
# max_z <- max(z_scaled_data)
# 
# scaled_data <- -1 + 2 * (z_scaled_data - min_z) / (max_z - min_z)
# 
# print("Scaled Data (-1 to 1):", scaled_data)


# 
# # 2. Z-Score Normalization to a Custom Range [-1, 1]
# min_z <- min(z_scaled_data)
# max_z <- max(z_scaled_data)
# scaled_data <- -1 + 2 * (z_scaled_data - min_z) / (max_z - min_z)
# print(scaled_data)
# 

RNA=read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Annotation_RNA_TMM.csv", header = TRUE,check.names = F,row.names = 1)

Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylated_stan_norm.csv", header = TRUE,check.names = F)
colnames(Methylation)[1]="X"
aggregated_Methylation <- aggregate(Methylation[, -1], by = list(X = Methylation$X), FUN = mean)
aggregated_Methylation=aggregated_Methylation[-c(1:10),]
rownames(aggregated_Methylation)=aggregated_Methylation[,1]
aggregated_Methylation <- aggregated_Methylation %>% select(-X)

aggregated_rna <- RNA %>%
  group_by(GeneSymbol) %>%
  summarise(across(.cols = starts_with("GSM"), .fns = mean, na.rm = TRUE))

aggregated_rna=as.data.frame(aggregated_rna)
aggregated_rna <- aggregated_rna[!is.na(aggregated_rna$GeneSymbol), ]
rownames(aggregated_rna)=aggregated_rna$GeneSymbol

aggregated_rna <- aggregated_rna %>% select(-GeneSymbol)
write.csv(aggregated_rna,"mapped_rna_TMM_normalized_scales.csv")

metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

# Rename columns in dataset if they exist in the metadata
colnames(aggregated_Methylation) <- rename_vector[colnames(aggregated_Methylation)]
colnames(aggregated_rna) <- rename_vector[colnames(aggregated_rna)]


aggregated_rna=aggregated_rna[,colnames(aggregated_Methylation)]

run_MOFA_model = function(omic1, omic2){
  #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
  #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
  
  #Creat_sample_Mofa_data_frame
  data <- make_example_data(
    n_views = 2,
    n_samples = 2,
    n_features = 2,
    n_factors = 2
  )[[1]]
  names(data)<-c("Transcriptomics","Methylatione")
  
  #Insert_omics_into_MOFAaccepted_data_frame
  RNA2 = as.data.frame(omic1)
  Methylatione = as.data.frame(omic2)
  #Convert_the_input_into_matrix
  omic1<-data.matrix(RNA2)
  omic2<-data.matrix(Methylatione)
  data[["Transcriptomics"]]<-omic1
  data[["Methylatione"]]<-omic2
  N = ncol(data[[1]])
  #MOFAobject <- create_mofa(data, groups=groups)
  MOFAobject <- create_mofa(data)
  #3.3 Visualise the structure of the data
  plot_data_overview(MOFAobject)
  # Define options
  # 4.1Define data options
  data_opts <- get_default_data_options(MOFAobject)
  data_opts$scale_views = T
  # data_opts$scale_views = T
  #4.2 Define model options
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- 5
  #4.3 Define train options
  train_opts <- get_default_training_options(MOFAobject)
  #5 Build and train the MOFA object
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile = "model16.hdf5"
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
  
  model <- load_model(outfile)
  return(model)
}


library(reticulate)
# In terminal
#pip install scikit-learn
#set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
# pip install mofapy2


metadata_new=metadata[1:18,]
rownames(aggregated_Methylation)=paste0(rownames(aggregated_Methylation),"_Methyl")
run_MOFA_model(aggregated_rna, aggregated_Methylation)
model_scaled = load_model("model16.hdf5",remove_inactive_factors = FALSE)
plot_data_overview(model_scaled)
sample_metadata <- model_scaled@samples_metadata
sample_metadata$run = metadata_new$Age

library(tidyr)
library(dplyr)
metadata_new <- metadata_new %>%
  separate(Title, into = c("Patient", "sample.type"), sep = ",")


sample_metadata$sample.type = metadata_new$sample.type
sample_metadata$Age = metadata_new$Age
sample_metadata$gen= metadata_new$Gender
samples_metadata(model_scaled) <- sample_metadata

plot_data_overview(model_scaled)
plot_variance_explained(model_scaled, x="view", y="factor", factors = c(1:5))

plot_variance_explained(model_scaled, x="group", y="factor", plot_total = T)[[2]]


plot_factor(model_scaled, 
            factor = 1:5,
            color_by = "sample.type")

plot_factor(model_scaled, 
            factor = 1:2,
            color_by = "gen")