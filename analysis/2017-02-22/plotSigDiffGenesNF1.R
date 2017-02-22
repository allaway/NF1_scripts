library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(ggplot2)
synapseLogin()

this.file<-

patients<-read.table("../../bin/TCGA_NF1_patient_status.txt")
patients<-filter(patients, mutation==1)
mutpts<-as.character(patients$samples)

##get muts
mutation.files<-list.files("../../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../../data/mutations_by_cancer/",x, sep=""), header = TRUE))
cancer.type<-sub("mutations.txt", "", mutation.files)
names(mutations) <- cancer.type


library(wesanderson)
library(viridis)

##get filenames for 
files <-
  synQuery("SELECT * from file WHERE parentId=='syn8303919'")$file.id

i <- "syn8303923"

for(i in files){
    data<-read.table(synGet(i)@filePath, sep = "\t", header = TRUE)
    cancer<-synGet(i)@fileHandle$fileName
    cancer<-sub("_NF1_comut_TCGA.txt", "", cancer)
    bar<-filter(data, p_value <= 0.01)
    list<-bar$gene
    dat<-filter(mutations[[cancer]], gene %in% list)
    dat<-select(dat, gene, one_of(mutpts))
    dat$sums<-rowSums(dat[,-1])
    dat<-select(dat, gene, sums)
    dat$gene <- factor(dat$gene, levels = dat$gene[order(desc(dat$sums))])
    virid <- rev(viridis(n = nrow(dat)))
    ggplot(data = dat, aes(x=gene, y=sums)) +
      geom_bar(stat = "identity", color = "darkgray", fill = virid) +
      theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
      labs(title = paste("NF1 co-mutated genes in ", cancer, " (corr. p<0.01)", sep = ""), x = "Gene", y = "Comutation Frequency") +
      theme(plot.title=element_text(hjust=0.5))
    ggsave(paste("sig_NF1_co_mutant_genes_in_",cancer,"_TCGA.png", sep = ""), plot = last_plot(), width = 7, height = 5.5) 
    synStore()
    }

