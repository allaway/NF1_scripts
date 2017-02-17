library(data.table)
library(synapseClient)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(mHG)
library(parallel)
cpu<-detectCores()

##read in mutation files previously processed from Firehose data 
mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE))
cancers <- sub("mutations.txt", "", mutation.files)
names(mutations) <- cancers

pancan<-mutations %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene"), .)
pancan<-pancan[, -grep(".y", ignore.case = FALSE, colnames(pancan))]
pancan[is.na(pancan)] <- 0
pancan$sums<-rowSums(pancan[,-1])
pancan.muts<-select(pancan, gene, sums)
pancan.muts<-arrange(pancan.muts, desc(sums))

mHGs<-mclapply(cancers, function(x){
  print(x)
  foo<-mutations[[x]]
  foo$sums<-rowSums(foo[,-1])
  bar<-select(foo, gene, sums)
  bar<-arrange(bar, desc(sums))
  N<-nrow(pancan.muts)+nrow(bar)
  B<-nrow(bar)
  lambdas<-numeric(N)
  lambdas[sample(N, B)] <- 1
  p<-mHG.statistic.calc(lambdas)@mHG
  mHG.pval.calc(p, N, B)
}, mc.cores=cpu)

mnames(mHGs) <- cancers
Hy