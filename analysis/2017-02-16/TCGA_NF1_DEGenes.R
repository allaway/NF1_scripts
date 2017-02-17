library(synapseClient)
library(plyr)
library(dplyr)
library(limma)
library(edgeR)


this.file <- "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-02-16/TCGA_NF1_DEGenes.R"

tcga<-getDisExpressionData()

nf1.stat.df<-read.table(synGet("syn8265248")@filePath, sep = "\t", header = TRUE)
nf1.stat.l<-dlply(nf1.stat.df, .var = "cancer")

for(i in names(nf1.stat.l)){
  try({
    print(i)
    nf1mut<-filter(nf1.stat.l[[i]], mutation == "1")
    nf1mut<-as.character(nf1mut[,"samples"])
    nf1wt<-filter(nf1.stat.l[[i]], mutation == "0")
    nf1wt<-as.character(nf1wt[,"samples"])
    
    nf1mut2<-nf1mut[nf1mut %in% colnames(tcga)]
    mutpts<-dplyr::select(tcga, one_of(nf1mut2))
    
    nf1wt2<-nf1wt[nf1wt %in% colnames(tcga)]
    wtpts<-dplyr::select(tcga, one_of(nf1wt2))
    
    pts<-round(cbind(mutpts, wtpts))
    
    annot<-dplyr::filter(nf1.stat.l[[i]], samples %in% colnames(pts) )
    annot$mutation<-as.factor(sub(1, "mutant", annot$mutation))
    annot$mutation<-as.factor(sub(0, "wt", annot$mutation))
    
    rownames(annot) <- annot$samples
    annot<-dplyr::select(annot, mutation)
    design<-model.matrix(~annot$mutation)
    pts <- pts[, rownames(annot)]
    
    pts <- pts[complete.cases(pts),]
    
    dge<-DGEList(counts=pts)
    v<-voom(dge, design, plot = TRUE)
    fit <- lmFit(v, design)
    fit <- eBayes(fit)
    res<-topTable(fit,sort="none",n=Inf)
    file<-write.table(res, file = paste("TCGA_",i,"_NF1_DEgenes.csv", sep = ""), sep = ",")
    synStore(File(paste("TCGA_",i,"_NF1_DEgenes.csv", sep = ""), parentId = "syn8267685"), executed = this.file)
  })
}



