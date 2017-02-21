library(data.table)
library(dplyr)
library(reshape2)
library(synapseClient)

cancers<-list.dirs("../data")
cancers<-sub("../data/tcga_mafs/gdac.broadinstitute.org_", "", cancers)
cancers<-sub(".Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0", "", cancers)
cancers<-sub(".Mutation_Packager_Oncotated_Raw_Calls.Level_3.2016012800.0.0", "", cancers)
cancers<-cancers[4:length(cancers)]
print(cancers)

mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE))

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
    g<-mutations[[i]][j,1]
    bar<-filter(mutations[[i]], gene == paste("\'",g,"\'", sep=""))
    bar2<-filter(mutations[[i]], gene == "NF1")
    bar3<-rbind(bar, bar2)
    bar3<-t(select(bar3, -sums))
    colnames(bar3) <- bar3[1,]
    bar3 <- as.data.frame(bar3[-1,])

    ptswNF1only<-filter(bar3, "NF1"==1)
    ptswNF1only<-nrow(filter_(ptswNF1only, paste("\'",g,"\'","==0", sep="")))

    ptswGonly<-filter(bar3, "NF1"==0)
    ptswGonly<-nrow(filter_(ptswGonly, paste("\'",g,"\'","==1", sep="")))
    
    ptswGandNF1<-filter(bar3, "NF1"==1)
    ptswGandNF1<-nrow(filter_(ptswGandNF1, paste("\'",g,"\'","==1", sep="")))
    
    nomut<-filter(bar3, "NF1"==0)
    nomut<-nrow(filter_(nomut, paste("\'",g,"\'","==0", sep="")))
    
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

cancer.type<-sub("mutations.txt", "", mutation.files)
names(comuts) <- cancer.type

comuts2<-lapply(names(comuts), function(k){
  if(nrow(comuts[[k]])>1){
  foo<-as.data.frame(comuts[[k]])
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
