

library(MOFA2)

RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Filtered_RNA_standard_normalized.csv" ,header = TRUE,check.names = F,row.names = 1)
rownames(RNA)=RNA$GeneSymbol
RNA <- RNA %>% select(-GeneSymbol)

# 
# Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylated_stan_norm.csv", header = TRUE,check.names = F)
# colnames(Methylation)[1]="X"
# aggregated_Methylation <- aggregate(Methylation[, -1], by = list(X = Methylation$X), FUN = mean)
# aggregated_Methylation=aggregated_Methylation[-c(1:10),]
# rownames(aggregated_Methylation)=aggregated_Methylation[,1]
# aggregated_Methylation <- aggregated_Methylation %>% select(-X)


metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

# Rename columns in dataset if they exist in the metadata
# colnames(aggregated_Methylation) <- rename_vector[colnames(aggregated_Methylation)]
colnames(RNA) <- rename_vector[colnames(RNA)]


# RNA=RNA[,colnames(aggregated_Methylation)]

run_MOFA_model = function(omic1){
  #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row

  #Creat_sample_Mofa_data_frame
  data <- make_example_data(
    n_views = 1,
    n_samples = 2,
    n_features = 2,
    n_factors = 2
  )[[1]]
  names(data)<-c("Transcriptomics")
  
  #Insert_omics_into_MOFAaccepted_data_frame
  RNA2 = as.data.frame(omic1)
  #Convert_the_input_into_matrix
  omic1<-data.matrix(RNA2)
  data[["Transcriptomics"]]<-omic1
  N = ncol(data[[1]])
  #MOFAobject <- create_mofa(data, groups=groups)
  MOFAobject <- create_mofa(data)
  #3.3 Visualise the structure of the data
  plot_data_overview(MOFAobject)
  # Define options
  # 4.1Define data options
  data_opts <- get_default_data_options(MOFAobject)
  # data_opts$scale_views = T
  #4.2 Define model options
  model_opts <- get_default_model_options(MOFAobject)
  
  #4.3 Define train options
  train_opts <- get_default_training_options(MOFAobject)
  #5 Build and train the MOFA object
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile = "model11.hdf5"
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
run_MOFA_model(RNA)
model = load_model("model11.hdf5")
plot_data_overview(model)
sample_metadata <- model@samples_metadata
sample_metadata$run = metadata_new$Age

library(tidyr)
library(dplyr)
metadata_new <- metadata_new %>%
  separate(Title, into = c("Patient", "sample.type"), sep = ",")


sample_metadata$sample.type = metadata_new$sample.type
sample_metadata$Age = metadata_new$Age
sample_metadata$gen= metadata_new$Gender
samples_metadata(model) <- sample_metadata

plot_data_overview(model)
plot_variance_explained(model, x="view", y="factor", factors = c(1:5))

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]


plot_factor(model, 
            factor = 1:5,
            color_by = "sample.type")

plot_factor(model, 
            factor = 1:5,
            color_by = "gen")

plot_data_heatmap(model,
                  view = "Transcriptomics",         # view of interest
                  factor = 7,             # factor of interest
                  features = 25,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the pheatmap function
                  cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
                  
                  show_rownames = TRUE, show_colnames = FALSE)

plot_data_heatmap(model,
                  view = "Methylatione",         # view of interest
                  factor = 4,             # factor of interest
                  features = 28,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the pheatmap function
                  cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
                  
                  show_rownames = FALSE, show_colnames = FALSE
)



