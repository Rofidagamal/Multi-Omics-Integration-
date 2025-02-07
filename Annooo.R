# df=read.csv("C:\\Users\\Electronica Care\\Downloads\\Telegram Desktop\\Annotated_Methylation_Data_DEMs.csv")
# 
# aggregated_data <- mapped_data %>%
#   group_by(UCSC_RefGene_Name) %>%
#   summarise(across(.cols = c( "lfc.diff" , "t.pval" , "t.pval.adj"), .fns = mean, na.rm = TRUE))
# 
# write.csv(aggregated_data,"DEMs.csv")

