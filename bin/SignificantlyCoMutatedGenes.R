library(data.table)
library(dplyr)
library(reshape2)
library(synapseClient)

mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE))

cancer.type<-sub("mutations.txt", "", mutation.files)
names(mutations) <- cancer.type

comuts<-lapply(1:35, function(i){
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
  foo<-lapply(mutations[[i]]$gene, function(j) {
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
    
  })
    
  #comuts<-unname(comuts)
  names(foo) <- mutations[[i]]$gene
  return(foo)
  
  } else {
    foo <- "no NF1 mutations"
    return (foo)
  }
  
})

comuts2<-lapply(comuts, function(k){
  if(length(k)>1){
  foo<-as.data.frame(k)
  bar<-t(foo)
  colnames(bar) <- c("p_value", "estimate")
  bar<-as.data.frame(bar)
  bar$BH<-p.adjust(bar[,1], method = "BH")
  return(bar)
  } else {
    bar <- "no nf1 mutations"
    return(bar)
  }
})

names(comuts2) <- cancer.type

pancan<-mutations %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene"), .)

pancan<-pancan[, -grep(".y", ignore.case = FALSE, colnames(pancan))]
pancan[is.na(pancan)] <- 0
