# Load required library
library(dplyr)
# # Read the CSV file (assuming it's saved as "data.csv")
# metadata2 <- metadata
# # Extract Patient ID from Title
# df <- df %>%
#   mutate(Patient_ID = gsub(",.*", "", Title)) %>%  # Extract "Patient X" part
#   group_by(Patient_ID) %>%
#   mutate(Accession_Mapped = first(Accession))  # Assign the same accession to lesion and NAWM
#
# # View the transformed data
# print(df)

library(MOFA2)

RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\RNA_normalized_counts_filtered.csv", header = TRUE, row.names = 1)
Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\dmrs_genes_count_data (1).csv", header = TRUE, row.names = 1)
metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

# Rename columns in dataset if they exist in the metadata
colnames(Methylation) <- rename_vector[colnames(Methylation)]
colnames(RNA) <- rename_vector[colnames(RNA)]

RNA=RNA[,colnames(Methylation)]

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

  outfile = "model.hdf5"
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
run_MOFA_model(RNA, Methylation)
model = load_model("model.hdf5")
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


p <- plot_factor(model, 
                 factors = c(1,2,3,4,5),
                 color_by = "sample.type",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
library(GGally)

p <- p + 
  scale_color_manual(values=c("lesion"="red","NAWM"="green")) +
  scale_fill_manual(values=c("lesion"="red","NAWM"="green"))

print(p)






# Load required library
library(dplyr)


# # Read the CSV file (assuming it's saved as "data.csv")
# metadata2 <- metadata
# # Extract Patient ID from Title
# df <- df %>%
#   mutate(Patient_ID = gsub(",.*", "", Title)) %>%  # Extract "Patient X" part
#   group_by(Patient_ID) %>%
#   mutate(Accession_Mapped = first(Accession))  # Assign the same accession to lesion and NAWM
# 
# # View the transformed data
# print(df)


#' 
#' 
#' 
#' 
#' 
#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_raw_counts_GRCh38.p13_NCBI (1).tsv", header = TRUE, row.names = 1,sep='\t')
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Methylation_20000.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model2.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model2.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:5))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 





#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\normalized_RNA.csv",header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Methylation_20000.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model3.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model3.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:5))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 
#' 







#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\GSE224377_norm_counts_FPKM_GRCh38.p13_NCBI.tsv",sep='\t',header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Methylation_20000.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   # data_opts$scale_views = T
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model4.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model4.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:5))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:4,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:4,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 
#' 





#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\filtered_dataset_RNA_less than 0.5.csv"
#' ,header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Methylation_20000.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   # data_opts$scale_views = T
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model5.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model5.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:5))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 


#' 
#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\RNA_filtered_means.csv" ,header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylation_filtered_means.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   # data_opts$scale_views = T
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model6.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model6.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:5))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 
#' 



#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\RNA_filtered_means.csv" ,header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylation_filtered_means.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   # data_opts$scale_views = T
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model7.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model7.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:5))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:5,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 



#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\RNA_filtered_means.csv" ,header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylation_filtered_means.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   # data_opts$scale_views = T
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   model_opts$num_factors <- 7
#' 
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model8.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model8.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:7))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:7,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:7,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 
#' 



#' 
#' library(MOFA2)
#' 
#' RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\RNA_filtered_means.csv" ,header = TRUE, row.names = 1,check.names = F)
#' Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylation_filtered_means.csv", header = TRUE, row.names = 1)
#' metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")
#' 
#' rename_vector <- setNames(metadata$X, metadata$Accession)
#' 
#' # Rename columns in dataset if they exist in the metadata
#' colnames(Methylation) <- rename_vector[colnames(Methylation)]
#' colnames(RNA) <- rename_vector[colnames(RNA)]
#' 
#' RNA=RNA[,colnames(Methylation)]
#' 
#' run_MOFA_model = function(omic1, omic2){
#'   #' @param omic1 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   #' @param omic2 data_frame, include_omic1_sample_in_coloumn_and_features_in_row
#'   
#'   #Creat_sample_Mofa_data_frame
#'   data <- make_example_data(
#'     n_views = 2,
#'     n_samples = 2,
#'     n_features = 2,
#'     n_factors = 2
#'   )[[1]]
#'   names(data)<-c("Transcriptomics","Methylatione")
#'   
#'   #Insert_omics_into_MOFAaccepted_data_frame
#'   RNA2 = as.data.frame(omic1)
#'   Methylatione = as.data.frame(omic2)
#'   #Convert_the_input_into_matrix
#'   omic1<-data.matrix(RNA2)
#'   omic2<-data.matrix(Methylatione)
#'   data[["Transcriptomics"]]<-omic1
#'   data[["Methylatione"]]<-omic2
#'   N = ncol(data[[1]])
#'   #MOFAobject <- create_mofa(data, groups=groups)
#'   MOFAobject <- create_mofa(data)
#'   #3.3 Visualise the structure of the data
#'   plot_data_overview(MOFAobject)
#'   # Define options
#'   # 4.1Define data options
#'   data_opts <- get_default_data_options(MOFAobject)
#'   # data_opts$scale_views = T
#'   #4.2 Define model options
#'   model_opts <- get_default_model_options(MOFAobject)
#'   model_opts$num_factors <- 11
#'   
#'   #4.3 Define train options
#'   train_opts <- get_default_training_options(MOFAobject)
#'   #5 Build and train the MOFA object
#'   MOFAobject <- prepare_mofa(
#'     object = MOFAobject,
#'     data_options = data_opts,
#'     model_options = model_opts,
#'     training_options = train_opts
#'   )
#'   
#'   outfile = "model9.hdf5"
#'   MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = F)
#'   
#'   model <- load_model(outfile)
#'   return(model)
#' }
#' 
#' 
#' library(reticulate)
#' # In terminal
#' #pip install scikit-learn
#' #set SKLEARN_ALLOW_DEPRECATED_SKLEARN_PACKAGE_INSTALL=True
#' # pip install mofapy2
#' 
#' 
#' metadata_new=metadata[1:18,]
#' run_MOFA_model(RNA, Methylation)
#' model = load_model("model9.hdf5")
#' plot_data_overview(model)
#' sample_metadata <- model@samples_metadata
#' sample_metadata$run = metadata_new$Age
#' 
#' library(tidyr)
#' library(dplyr)
#' metadata_new <- metadata_new %>%
#'   separate(Title, into = c("Patient", "sample.type"), sep = ",")
#' 
#' 
#' sample_metadata$sample.type = metadata_new$sample.type
#' sample_metadata$Age = metadata_new$Age
#' sample_metadata$gen= metadata_new$Gender
#' samples_metadata(model) <- sample_metadata
#' 
#' plot_data_overview(model)
#' plot_variance_explained(model, x="view", y="factor", factors = c(1:11))
#' 
#' plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
#' 
#' 
#' plot_factor(model, 
#'             factor = 1:7,
#'             color_by = "sample.type")
#' 
#' plot_factor(model, 
#'             factor = 1:7,
#'             color_by = "gen")
#' 
#' plot_data_heatmap(model,
#'                   view = "Transcriptomics",         # view of interest
#'                   factor = 7,             # factor of interest
#'                   features = 25,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = FALSE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = TRUE, show_colnames = FALSE)
#' 
#' plot_data_heatmap(model,
#'                   view = "Methylatione",         # view of interest
#'                   factor = 4,             # factor of interest
#'                   features = 28,          # number of features to plot (they are selected by weight)
#'                   
#'                   # extra arguments that are passed to the pheatmap function
#'                   cluster_rows = TRUE, cluster_cols = TRUE,annotation_samples = "infection",
#'                   
#'                   show_rownames = FALSE, show_colnames = FALSE
#' )
#' 
#' 




library(MOFA2)

RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Filtered_RNA_standard_normalized.csv" ,header = TRUE,check.names = F,row.names = 1)
rownames(RNA)=RNA$GeneSymbol
RNA <- RNA %>% select(-GeneSymbol)


Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylated_stan_norm.csv", header = TRUE,check.names = F)
colnames(Methylation)[1]="X"
aggregated_Methylation <- aggregate(Methylation[, -1], by = list(X = Methylation$X), FUN = mean)
aggregated_Methylation=aggregated_Methylation[-c(1:10),]
rownames(aggregated_Methylation)=aggregated_Methylation[,1]
aggregated_Methylation <- aggregated_Methylation %>% select(-X)


metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

# Rename columns in dataset if they exist in the metadata
colnames(aggregated_Methylation) <- rename_vector[colnames(aggregated_Methylation)]
colnames(RNA) <- rename_vector[colnames(RNA)]


RNA=RNA[,colnames(aggregated_Methylation)]

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
  
  outfile = "model10.hdf5"
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
run_MOFA_model(RNA, aggregated_Methylation)
model_standard_norm = load_model("model10.hdf5")
plot_data_overview(model_standard_norm)
sample_metadata <- model_standard_norm@samples_metadata
sample_metadata$run = metadata_new$Age

library(tidyr)
library(dplyr)
metadata_new <- metadata_new %>%
  separate(Title, into = c("Patient", "sample.type"), sep = ",")


sample_metadata$sample.type = metadata_new$sample.type
sample_metadata$Age = metadata_new$Age
sample_metadata$gen= metadata_new$Gender
samples_metadata(model_standard_norm) <- sample_metadata

plot_data_overview(model_standard_norm)
plot_variance_explained(model_standard_norm, x="view", y="factor", factors = c(1:5))

plot_variance_explained(model_standard_norm, x="group", y="factor", plot_total = T)[[2]]


plot_factor(model_standard_norm, 
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



















RNA=as.matrix(RNA)
# 1. Standard Z-Scaling (Z-score normalization)
mu <- mean(RNA)
sigma <- sd(RNA)
z_scaled_data <- (RNA - mu) / sigma



# 2. Z-Scaling and then Min-Max Scaling to -1 to 1

# (Z-scaling is already done above, so we just reuse z_scaled_data)

min_z <- min(z_scaled_data)
max_z <- max(z_scaled_data)

scaled_data <- -1 + 2 * (z_scaled_data - min_z) / (max_z - min_z)

print("Scaled Data (-1 to 1):", scaled_data)




library(MOFA2)

RNA= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Filtered_RNA_standard_normalized.csv" ,header = TRUE,check.names = F,row.names = 1)
rownames(RNA)=RNA$GeneSymbol
RNA <- RNA %>% select(-GeneSymbol)


Methylation= read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Methylated_stan_norm.csv", header = TRUE,check.names = F)
colnames(Methylation)[1]="X"
aggregated_Methylation <- aggregate(Methylation[, -1], by = list(X = Methylation$X), FUN = mean)
aggregated_Methylation=aggregated_Methylation[-c(1:10),]
rownames(aggregated_Methylation)=aggregated_Methylation[,1]
aggregated_Methylation <- aggregated_Methylation %>% select(-X)


metadata = read.csv("C:\\Users\\Electronica Care\\Downloads\\Metadata_RNA_Methyl.csv")

rename_vector <- setNames(metadata$X, metadata$Accession)

# Rename columns in dataset if they exist in the metadata
colnames(aggregated_Methylation) <- rename_vector[colnames(aggregated_Methylation)]
colnames(RNA) <- rename_vector[colnames(RNA)]


RNA=RNA[,colnames(aggregated_Methylation)]

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
  
  outfile = "model13.hdf5"
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
run_MOFA_model(scaled_data, aggregated_Methylation)
model_scaled = load_model("model13.hdf5")
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

plot_data_overview(model)
plot_variance_explained(model_scaled, x="view", y="factor", factors = c(1:5))

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]


plot_factor(model, 
            factor = 1:5,
            color_by = "sample.type")

plot_factor(model, 
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



