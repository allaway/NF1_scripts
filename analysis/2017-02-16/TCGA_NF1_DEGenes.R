library(synapseClient)
library(plyr)
library(dplyr)
library(DESeq2)


tcga<-getDisExpressionData()

nf1.stat.df<-read.table(synGet("syn8265248")@filePath, sep = "\t", header = TRUE)
nf1.stat.l<-dlply(nf1.stat.df, .var = "cancer")

for(i in names(nf1.stat.l, function(i)){
  try({
    print(i)
    nf1mut<-filter(nf1.stat.l[[i]], mutation == "1")
    nf1mut<-as.character(nf1mut[,"samples"])
    nf1wt<-filter(nf1.stat.l[[i]], mutation == "0")
    nf1wt<-as.character(nf1wt[,"samples"])
    
    nf1mut2<-nf1mut[nf1mut %in% colnames(tcga)]
    mutpts<-select(tcga, one_of(nf1mut2))
    
    nf1wt2<-nf1wt[nf1wt %in% colnames(tcga)]
    wtpts<-select(tcga, one_of(nf1wt2))
    
    pts<-round(cbind(mutpts, wtpts))
    
    annot<-filter(nf1.stat.l[[i]], samples %in% colnames(pts) )
    annot$mutation<-as.factor(sub(1, "mutant", annot$mutation))
    annot$mutation<-as.factor(sub(0, "wt", annot$mutation))
    
    rownames(annot) <- annot$samples
    annot<-select(annot, mutation)
    pts <- pts[, rownames(annot)]
    
    pts <- pts[complete.cases(pts),]
    
    countData<-DESeqDataSetFromMatrix(countData = pts, colData = annot, design = ~ mutation)
    res <- DESeq(countData)
    file<-write.table(as.data.frame(res), file = paste("TCGA_",i,"_NF1_DEgenes.csv", sep = ""), sep = ",")
    synStore(File(file, parentId = ))
  })
}


