library(MOFA2)

Rna_counts <- read.csv("/Users/hageradel/Downloads/GSE224377_full_table_counts_QC_normalized_CSV.csv", header = T, row.names=1, check.names=F)
Methylation= read.csv("/Users/hageradel/Downloads/GSE224455_Methylated_intensities_CSV.csv", header = T, row.names=1, check.names=F)
metadata = read.csv("/Users/hageradel/Desktop/Data/Methyl_RNA.data - Sheet1 (1).csv", header = T, check.names=F)

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
  names(data)<-c("Transcriptomics","Methylation")
  
  #Insert_omics_into_MOFAaccepted_data_frame
  Rna_counts = as.data.frame(omic1)
  Methylation = as.data.frame(omic2)
  #Convert_the_input_into_matrix
  omic1<-data.matrix(Rna_counts)
  omic2<-data.matrix(Methylation)
  data[["Transcriptomics"]]<-omic1
  data[["Methylation"]]<-omic2
  N = ncol(data[[1]])
  #MOFAobject <- create_mofa(data, groups=groups)
  MOFAobject <- create_mofa(data)
  #3.3 Visualise the structure of the data
  plot_data_overview(MOFAobject)
  # Define options
  # 4.1Define data options
  data_opts <- get_default_data_options(MOFAobject)
  #4.2 Define model options
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors=7
  #4.3 Define train options
  train_opts <- get_default_training_options(MOFAobject)
  #5 Build and train the MOFA object
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile = "/Users/hageradel/Downloads/SOLOmodel.hdf5"
  MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = T)
  
  model <- load_model(outfile)
  return(model)
}


library(reticulate)
# In terminal
#pip install scikit-learn
#set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
# pip install mofapy2


run_MOFA_model(Rna_counts, Methylation)
model = load_model("/Users/hageradel/Downloads/SOLOmodel.hdf5")
plot_data_overview(model)
sample_metadata <- model@samples_metadata
sample_metadata$run = metadata$Age
sample_metadata$infection = metadata$infection
sample_metadata$Age = metadata$Age_cat
sample_metadata$gen= metadata$gender
samples_metadata(model) <- sample_metadata

plot_data_overview(model)
plot_variance_explained(model, x="view", y="factor", factors = c(1:7))

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]


plot_factor(model, 
            factor = 1:5,
            color_by = "infection")

plot_data_heatmap(model,
                  view = "Transcriptomics",         # view of interest
                  factor = 7,             # factor of interest
                  features = 25,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
                  
                  show_rownames = TRUE, show_colnames = FALSE)

plot_data_heatmap(model,
                  view = "Micro",         # view of interest
                  factor = 4,             # factor of interest
                  features = 28,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
                  
                  show_rownames = FALSE, show_colnames = FALSE
)
