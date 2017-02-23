library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(ggplot2)
synapseLogin()

this.file<-"https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-02-22/plotSigCoMutantGenesNF1.R"

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

for(i in files){
    data<-read.table(synGet(i)@filePath, sep = "\t", header = TRUE)
    cancer<-synGet(i)@fileHandle$fileName
    cancer<-sub("_NF1_comut_TCGA.txt", "", cancer)
    bar<-filter(data, p_value <= 0.1)
    list<-bar$gene
    dat<-filter(mutations[[cancer]], gene %in% list)
    dat<-select(dat, gene, one_of(mutpts))
    
    if(NCOL(dat[,-1])>1 & NROW(dat[,-1])>1){
    dat$sums<-rowSums(dat[,-1])
    dat<-select(dat, gene, sums)
    dat$gene <- factor(dat$gene, levels = dat$gene[order(desc(dat$sums))])
    
    if(length(dat$gene)>30){
      dat<-arrange(dat, desc(sums))
      dat<-slice(dat, 1:30)
    }
    virid <- rev(viridis(n = nrow(dat)))
    ggplot(data = dat, aes(x=gene, y=sums)) +
      geom_bar(stat = "identity", color = "darkgray", fill = virid) +
      theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
      labs(title = paste("NF1 co-mutated genes in ", cancer, " (corr. p<0.1, top 30)", sep = ""), x = "Gene", y = "Comutation Frequency") +
      theme(plot.title=element_text(hjust=0.5))
    ggsave(paste("sig_NF1_co_mutant_genes_in_",cancer,"_TCGA.png", sep = ""), plot = last_plot(), width = 7, height = 5.5) 
    synStore(File(paste("sig_NF1_co_mutant_genes_in_",cancer,"_TCGA.png", sep = ""),parentId = "syn8304189"), used = i, executed = this.file) 
    } else {
      print("no sig comuts") 
    }
    
    }

