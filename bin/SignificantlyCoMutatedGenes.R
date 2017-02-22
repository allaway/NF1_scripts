library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(ggplot2)
library(parallel)

synapseLogin()

this.file= "https://raw.githubusercontent.com/allaway/NF1_scripts/master/bin/SignificantlyCoMutatedGenes.R"

mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE))

cancer.type<-sub("mutations.txt", "", mutation.files)
names(mutations) <- cancer.type

comuts<-mclapply(1:35, function(i){
  ##get samples with NF1 mutation
  nf1<-dplyr::filter(mutations[[i]], gene=="NF1")
  nf1<-nf1[,-1]
  nf1<-as.data.frame(t(nf1))
  nf1$samples<-rownames(nf1)
  nf1<-dplyr::filter(nf1, nf1[,1]==1)
  
  if(nrow(nf1)>0){
    cancer.type<-sub("mutations.txt", "", mutation.files[[i]])
    print(cancer.type)
    
    ##get all mutation sums 
    mutations[[i]]$sums<-rowSums(mutations[[i]][,-1])
    
    ##hypergeometric test for each gene
    foo<-sapply(mutations[[i]]$gene, function(j) {
      bar<-filter(mutations[[i]], gene == `j`)
      bar2<-filter(mutations[[i]], gene == "NF1")
      bar3<-t(rbind(bar, bar2))
      colnames(bar3) <- bar3[1,]
      bar3 <- as.data.frame(bar3[-1,])
      
      if(j=="NF1"){
        j="NF1.test"
        colnames(bar3)<-c("NF1", "NF1.test")
      }
      
      ptswNF1only<-filter(bar3, NF1==1)
      ptswNF1only<-nrow(filter_(ptswNF1only, paste("`",j,"`","==0", sep="")))
      
      ptswGonly<-filter(bar3, NF1==0)
      ptswGonly<-nrow(filter_(ptswGonly, paste("`",j,"`","==1", sep="")))
      
      ptswGandNF1<-filter(bar3, NF1==1)
      ptswGandNF1<-nrow(filter_(ptswGandNF1, paste("`",j,"`","==1", sep="")))
      
      nomut<-filter(bar3, NF1==0)
      nomut<-nrow(filter_(nomut, paste("`",j,"`","==0", sep="")))
      
      bar<-c(-1,-1)
      f<-fisher.test(matrix(c(ptswGandNF1,ptswNF1only,ptswGonly,nomut), ncol = 2))
      try(bar<-c(f$p.value, f$estimate))
      
      bar<-unname(bar) 
      
    }, simplify = FALSE, USE.NAMES = TRUE)
    
    #comuts<-unname(comuts)
    names(foo) <- mutations[[i]]$gene
    return(foo)
    
  } else {
    foo <- "no NF1 mutations"
    return (foo)
  }
  
}, mc.cores = detectCores())

comuts2<-mclapply(comuts, function(k) ldply(k), mc.cores = detectCores())

names(comuts2) <- cancer.type

comuts3<-lapply(comuts2, function(m){
  if(nrow(m)>1){
    bar<-m
    colnames(bar) <- c("gene","p_value","estimate")
    return(bar)
    
  } else {
    bar <- "no nf1 mutations"
    return(bar)
  }
})

for(q in names(comuts3)){
  print(q)
  if(length(comuts3[[q]])>1){
    write.table(comuts3[[q]], file = paste("../data/TCGA_Comuts/",q,"_NF1_comut_TCGA.txt",sep=""), sep = "\t")
    synStore(File(paste("../data/TCGA_Comuts/",q,"_NF1_comut_TCGA.txt",sep=""), parentId="syn8299578"), executed = this.file)
    foo<-filter(comuts3[[q]], "p_value" <= 0.01)
    write.table(foo, file = paste("../data/TCGA_Comuts/",q,"_NF1_comut_TCGA_010.txt",sep=""), sep = "\t")
    synStore(File(paste("../data/TCGA_Comuts/",q,"_NF1_comut_TCGA_010.txt",sep=""), parentId="syn8299578"), executed = this.file)         
    bar<-filter(comuts3[[q]], "p_value" <= 0.05)     
    write.table(bar, file = paste(./"data/TCGA_Comuts/",q,"_NF1_comut_TCGA_005.txt",sep=""), sep = "\t")
    synStore(File(paste("../data/TCGA_Comuts/",q,"_NF1_comut_TCGA_005.txt",sep=""), parentId="syn8299578"), executed = this.file)  
  }else{
    print("no NF1 mutation")
  }
}