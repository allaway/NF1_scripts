library(plyr)
library(dplyr)
library(synapseClient)
library(parallel)

synapseLogin()

this.file= "https://raw.githubusercontent.com/allaway/NF1_scripts/master/bin/SignificantlyCoMutatedGenesCCLE.R"

mutations<-read.table(synGet("syn7466552")@filePath)
gene<-rownames(mutations)
mutations<-cbind(gene, mutations)

nf1<-dplyr::filter(mutations, gene=="NF1")
nf1<-nf1[,-1]
nf1<-as.data.frame(t(nf1))
nf1$samples<-rownames(nf1)
nf1<-dplyr::filter(nf1, nf1[,1]==1)
  
##get all mutation sums 
mutations$sums<-rowSums(mutations[,-1])
    
##hypergeometric test for each gene
foo<-sapply(mutations$gene, function(j) {
  bar<-filter(mutations, gene == `j`)
  bar2<-filter(mutations, gene == "NF1")
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

names(foo) <- gene

comuts2 <- ldply(foo)
colnames(comuts2) <- c("gene","p_value","estimate")

write.table(comuts2, file = paste("../data/CCLE_NF1_comut.txt",sep=""), sep = "\t")
synStore(File(paste("../data/CCLE_NF1_comut.txt",sep=""), parentId="syn8495175"), executed = this.file, used = c("syn7466552"))  
