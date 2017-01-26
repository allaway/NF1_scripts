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

readInMafs<-function(cancerType){
  mani<-fread(file=paste("../data/tcga_mafs/gdac.broadinstitute.org_",cancerType,".Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/MANIFEST.txt", sep=""), header = FALSE)
  filenames<-mani$V2
  mut.df <- lapply(filenames, function(x) {
    dat <- try(fread(file = paste("../data/tcga_mafs/gdac.broadinstitute.org_",cancerType,".Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/",x, sep=""), 
                    skip = 3, header = TRUE, data.table = FALSE), silent = TRUE)
    dat <- dplyr::filter(dat, Variant_Classification != "Silent")
    dat <- dplyr::filter(dat, Variant_Classification !=  "Intron")
    mut.genes <- dplyr::select(dat, Hugo_Symbol)
    patient <- x
    patient <- sub(".hg19.oncotator.hugo_entrez_remapped.maf.txt", "", patient)
    mut.genes <- dplyr::distinct(mut.genes)
    obs<-rep.int(1, nrow(mut.genes))
    obs<-data.frame(obs)
    obs <- cbind(mut.genes, obs)
    colnames(obs) <- c("Hugo_Genes", patient)
    return(obs)
  })
  mut.df<-plyr::rbind.fill(mut.df)
  uniquegenes <- unique(mut.df$Hugo_Genes)
  mut.df[is.na(mut.df)] <- 0
  
 muts<-lapply(uniquegenes, function(i) {
    res<-dplyr::filter(mut.df, Hugo_Genes==i)
    res<-res[,-1]
    res<-colSums(res[,,drop=FALSE], na.rm = TRUE)
    res<-append(list("gene"=i), res)
    })
 
 muts<-rbindlist(muts, use.names=TRUE)
 write.table(muts, file = paste("../data/mutations_by_cancer/",cancerType,"mutations.txt", sep=""), sep = "\t")
}

for(m in cancers) {
  readInMafs(m)
  print(m)
}

mutation.files<-list.files("../data/mutations_by_cancer/")
mutations<-lapply(mutation.files, function(x) read.table(file = paste("../data/mutations_by_cancer/",x, sep=""), header = TRUE))

library(ggplot2)
library(wesanderson)
wes<-wes_palette("Darjeeling", 11, type = "continuous")

for(i in 1:35) {
  nf1<-dplyr::filter(mutations[[i]], gene=="NF2")
  nf1<-nf1[,-1]
  nf1<-as.data.frame(t(nf1))
  nf1$samples<-rownames(nf1)
  nf1<-dplyr::filter(nf1, nf1[,1]==1)
  cancer.type<-sub("mutations.txt", "", mutation.files[[i]])
  try(samples.with.nf1.mut<-dplyr::select(mutations[[i]], gene, one_of(nf1$samples)))
  try(samples.with.nf1.mut<-dplyr::filter(samples.with.nf1.mut, gene!="Unknown"))
  try(samples.with.nf1.mut$sums<-rowSums(samples.with.nf1.mut[,-1]))
  try(samples.with.nf1.mut<-samples.with.nf1.mut[order(-samples.with.nf1.mut$sums),])
  try(samples.with.nf1.mut$gene <- factor(samples.with.nf1.mut$gene, levels = samples.with.nf1.mut$gene))
  try(plot<-ggplot(samples.with.nf1.mut[1:11,], aes(x=samples.with.nf1.mut$gene[1:11], y=samples.with.nf1.mut$sums[1:11], fill=samples.with.nf1.mut$gene[1:11]))
      +geom_bar(stat="identity")
      +scale_fill_manual(values=wes)
      +labs(title=cancer.type, x = "Gene", y = "# co-mutated")
      +theme(axis.text.x = element_text(angle = 60, hjust = 1))
      +guides(fill = "none"))
  try(ggsave(paste(cancer.type,"-NF2.png", sep="")))
  try(write.table(file=paste("../data/NF2_only_samples/NF2only",cancer.type,"mutations.txt", sep=""), samples.with.nf1.mut, sep="\t"))
  #synStore()
}

pancan<-mutations %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="gene"), .)

pancan<-pancan[, -grep(".y", ignore.case = FALSE, colnames(pancan))]
pancan[is.na(pancan)] <- 0

###trying with other genes a-la carte
panwes<-wes_palette("Darjeeling", 30, type = "continuous")
nf1<-dplyr::filter(pancan, gene=="SUZ12")
nf1<-nf1[,-1]
nf1<-as.data.frame(t(nf1))
nf1$samples<-rownames(nf1)
nf1<-dplyr::filter(nf1, nf1[,1]==1)
cancer.type<-"pancan"
try(samples.with.nf1.mut<-dplyr::select(pancan, gene, one_of(nf1$samples)))
try(samples.with.nf1.mut<-dplyr::filter(samples.with.nf1.mut, gene!="Unknown"))
try(samples.with.nf1.mut$sums<-rowSums(samples.with.nf1.mut[,-1]))
try(samples.with.nf1.mut<-samples.with.nf1.mut[order(-samples.with.nf1.mut$sums),])
try(samples.with.nf1.mut$gene <- factor(samples.with.nf1.mut$gene, levels = samples.with.nf1.mut$gene))
try(plot<-ggplot(samples.with.nf1.mut[1:30,], aes(x=samples.with.nf1.mut$gene[1:30], y=samples.with.nf1.mut$sums[1:30], fill=samples.with.nf1.mut$gene[1:30]))
    +geom_bar(stat="identity")
    +scale_fill_manual(values=panwes)
    +labs(title=cancer.type, x = "Gene", y = "# co-mutated")
    +theme(axis.text.x = element_text(angle = 60, hjust = 1))
    +guides(fill = "none"))
try(ggsave("pancanSUZ12.png"))
