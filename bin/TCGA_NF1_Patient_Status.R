library(data.table)
library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)

mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE))

cancers<-sub("mutations.txt", "", mutation.files)

nf1.stat<-lapply(1:37, function(i){
  ##get samples with NF1 mutation
  nf1<-filter(mutations[[i]], gene=="NF1")
  nf1<-nf1[,-1]
  nf1<-as.data.frame(t(nf1))
  nf1$samples<-rownames(nf1)
  try(colnames(nf1) <- c("mutation", "samples"))
  cancer.type<-sub("mutations.txt", "", mutation.files[[i]])
  print(cancer.type)
  return(nf1)
})
names(nf1.stat) <- cancers
nf1.stat.df <- ldply(nf1.stat, .id = "cancer")
write.table(nf1.stat.df, file = "TCGA_NF1_patient_status.txt", sep = "\t")
synStore(File("TCGA_NF1_patient_status.txt", parentId = "syn7154899"))
