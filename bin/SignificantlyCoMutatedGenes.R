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

comuts<-lapply(1, function(i){
  ##get samples with NF1 mutation
  nf1<-dplyr::filter(mutations[[i]], gene=="NF1")
  nf1<-nf1[,-1]
  nf1<-as.data.frame(t(nf1))
  nf1$samples<-rownames(nf1)
  nf1<-dplyr::filter(nf1, nf1[,1]==1)
  cancer.type<-sub("mutations.txt", "", mutation.files[[i]])
  print(cancer.type)
  
  try({
  samples.with.nf1.mut<-dplyr::select(mutations[[i]], gene, one_of(nf1$samples))
  samples.with.nf1.mut<-dplyr::filter(samples.with.nf1.mut, gene!="Unknown")
  samples.with.nf1.mut$sums<-rowSums(samples.with.nf1.mut[,-1])
  #samples.with.nf1.mut<-samples.with.nf1.mut[order(-samples.with.nf1.mut$sums),]
  #samples.with.nf1.mut$gene <- factor(samples.with.nf1.mut$gene, levels = samples.with.nf1.mut$gene)
  })
  ##get all mutation sums 
  mutations[[i]]$sums<-rowSums(mutations[[i]][,-1])
  
  ##hypergeometric test for each gene
  foo<-lapply(mutations[[i]]$gene, function(j) {
    try({
    g<-mutations[[i]][j,1]
    bar<-filter(mutations[[i]], gene == g)
    bar2<-filter(mutations[[i]], gene == "NF1")
    bar3<-rbind(bar, bar2)
    bar3<-t(select(bar3, -sums))
    colnames(bar3) <- bar3[1,]
    bar3 <- as.data.frame(bar3[-1,])

    ptswNF1only<-filter(bar3, NF1==1)
    ptswNF1only<-nrow(filter_(ptswNF1only, paste(g,"==0")))

    ptswJonly<-filter(bar3, NF1==0)
    ptswJonly<-nrow(filter_(ptswJonly, paste(g,"==1")))
    
    ptswJandNF1<-filter(bar3, NF1==1)
    ptswJandNF1<-nrow(filter_(ptswJandNF1, paste(g,"==1")))
    
    nomut<-filter(bar3, NF1==0)
    nomut<-nrow(filter_(nomut, paste(g,"==0")))
    
    f<-fisher.test((matrix(c(ptswJandNF1,ptswNF1only,ptswJonly,nomut), ncol = 2)))
    
    bar<-c(f$p.value, f$estimate)
    bar<-unname(bar)
    if(is.numeric(bar)){
    return(bar)
    } else {
    bar<-c(-1,-1)
    }
    
  } , silent = TRUE) 
  })
  
  #comuts<-unname(comuts)
  names(foo) <- mutations[[i]]$gene
  return(foo)
  #ldply(foo)
  
})

cancer.type<-sub("mutations.txt", "", mutation.files)
names(muts) <- cancer.type

pancan<-mutations %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene"), .)

pancan<-pancan[, -grep(".y", ignore.case = FALSE, colnames(pancan))]
pancan[is.na(pancan)] <- 0
