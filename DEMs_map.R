DEMs=read.csv("../DEMs_final .csv")
DEMs_unique <- DEMs[!duplicated(DEMs$UCSC_RefGene_Name, fromLast = TRUE), ]
write.csv(DEMs_unique,"DEMs_Maps_unique.csv")

