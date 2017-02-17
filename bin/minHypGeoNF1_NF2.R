library(data.table)
library(dplyr)
library(reshape2)
library(synapseClient)
library(parallel)

mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-mclapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE), mc.cores = detectCores())
cancers <- sub("mutations.txt", "", mutation.files)
names(mutations) <- cancers

library(ggplot2)
library(wesanderson)
library(mHG)

for(i in cancers) {
  print(i)
  nf1<-dplyr::filter(mutations[[i]], gene=="NF1")
  nf1<-nf1[,-1]
  nf1<-as.data.frame(t(nf1))
  nf1$samples<-rownames(nf1)
  nf1<-dplyr::filter(nf1, nf1[,1]==1)
  cancer.type<-i
  try({
    samples.with.nf1.mut<-dplyr::select(mutations[[i]], gene, one_of(nf1$samples))
    samples.with.nf1.mut<-dplyr::filter(samples.with.nf1.mut, gene!="Unknown")
    samples.with.nf1.mut$sums<-rowSums(samples.with.nf1.mut[,-1])
    samples.with.nf1.mut<-samples.with.nf1.mut[order(-samples.with.nf1.mut$sums),]
    samples.with.nf1.mut$gene <- factor(samples.with.nf1.mut$gene, levels = samples.with.nf1.mut$gene)
    N<-nrow(mutations[[i]])
    B<-nrow(filter(samples.with.nf1.mut, sums > 0))
    lambdas<-numeric(N)
    lambdas[sample(N, B)] <- 1
    lambdas <- sort(lambdas, decreasing = TRUE)
    mhg<-mHG.test(lambdas)
    p<-mhg$p.value
    n<-mhg$n
    wes<-wes_palette("Darjeeling", n, type = "continuous")
    sigcomut <- samples.with.nf1.mut[1:n, ]
    plot<-ggplot(sigcomut, aes(x=sigcomut$gene, y=sigcomut$sums, fill=sigcomut$gene)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=wes) +
      labs(title=paste(cancer.type,"pval=",p, sep=" "), x = "Gene", y = "# co-mutated") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      guides(fill = "none")
    ggsave(paste(cancer.type,"-NF1.png", sep=""))
    write.table(file=paste("../data/NF1_only_samples/NF1cosig",cancer.type,"mutations.txt", sep=""), sigcomut, sep="\t")
  })
  #synStore()
}

pancan<-mutations %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene"), .)

pancan<-pancan[, -grep(".y", ignore.case = FALSE, colnames(pancan))]
pancan[is.na(pancan)] <- 0


