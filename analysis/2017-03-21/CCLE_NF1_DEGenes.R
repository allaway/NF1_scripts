library(plyr)
source("../../bin/cBioPortalData.R")
library(synapseClient)
library(dplyr)
library(limma)
library(edgeR)
synapseLogin()

this.file <- "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-03-21/CCLE_NF1_DEGenes.R"

##get expression data
ccle<-read.table(synGet("syn8369340")@filePath, sep = "\t", header = TRUE)

ccle.ov <- ccle[grep("*OVARY*", colnames(ccle))]
colnames(ccle.ov) <- sub("_.*","",colnames(ccle.ov))

ccle.brain <- ccle[grep("*CENTRAL_NERVOUS_SYSTEM*", colnames(ccle))]
colnames(ccle.brain) <- sub("_.*","",colnames(ccle.brain))

ccle<-select(ccle, -TT_OESOPHAGUS, -TT_THYROID)
colnames(ccle) <- sub("_.*","",colnames(ccle))

##get mutation data to ID NF1 mut samples
mutations<-read.table(synGet("syn7466552")@filePath)
gene<-rownames(mutations)
mutations<-cbind(gene, mutations)
nf1muts<-as.data.frame(t(filter(mutations, gene == "NF1")))
nf1muts$samps<-rownames(nf1muts)
nf1muts<-filter(nf1muts, V1 == "1")$samps

mut<-dplyr::select(ccle, one_of(nf1muts))
wt<-dplyr::select(ccle, -one_of(nf1muts))
    
pts<-round(cbind(mut,wt))

mutation<-c(rep("mutant", ncol(mut)), rep("wt", ncol(wt)))
names(mutation) <- c(colnames(mut), colnames(wt))
annot<-as.data.frame(mutation)

design<-model.matrix(~0+annot$mutation)
colnames(design)<-c("MU", "WT")
pts <- pts[, rownames(annot)]
pts <- pts[complete.cases(pts),]
    
fit <- lmFit(pts, design)
cont.matrix <- makeContrasts(MUvsWT=MU-WT, levels=design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
res<-topTable(fit, coef="MUvsWT", adjust="BH")

file<-write.table(res, file = paste("CCLE_NF1_DEgenes.csv", sep = ""), sep = ",")
synStore(File(paste("CCLE_NF1_DEgenes.csv", sep = ""), parentId = "syn8506635"), executed = this.file)


##just ovarian cancer
mut<-dplyr::select(ccle.ov, one_of(nf1muts))
wt<-dplyr::select(ccle.ov, -one_of(nf1muts))

pts<-round(cbind(mut,wt))

mutation<-c(rep("mutant", ncol(mut)), rep("wt", ncol(wt)))
names(mutation) <- c(colnames(mut), colnames(wt))
annot<-as.data.frame(mutation)

design<-model.matrix(~0+annot$mutation)
colnames(design)<-c("MU", "WT")
pts <- pts[, rownames(annot)]
pts <- pts[complete.cases(pts),]

fit <- lmFit(pts, design)
cont.matrix <- makeContrasts(MUvsWT=MU-WT, levels=design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
res<-topTable(fit, coef="MUvsWT", adjust="BH")

file<-write.table(res, file = paste("CCLE_OV_NF1_DEgenes.csv", sep = ""), sep = ",")
synStore(File(paste("CCLE_OV_NF1_DEgenes.csv", sep = ""), parentId = "syn8506635"), executed = this.file)


##just cns cancer
mut<-dplyr::select(ccle.brain, one_of(nf1muts))
wt<-dplyr::select(ccle.brain, -one_of(nf1muts))

pts<-round(cbind(mut,wt))

mutation<-c(rep("mutant", ncol(mut)), rep("wt", ncol(wt)))
names(mutation) <- c(colnames(mut), colnames(wt))
annot<-as.data.frame(mutation)

design<-model.matrix(~0+annot$mutation)
colnames(design)<-c("MU", "WT")
pts <- pts[, rownames(annot)]
pts <- pts[complete.cases(pts),]

fit <- lmFit(pts, design)
cont.matrix <- makeContrasts(MUvsWT=MU-WT, levels=design)
fit <- contrasts.fit(fit, cont.matrix)
fit <- eBayes(fit)
res<-topTable(fit, coef="MUvsWT", adjust="BH")

file<-write.table(res, file = paste("CCLE_CNS_NF1_DEgenes.csv", sep = ""), sep = ",")
synStore(File(paste("CCLE_CNS_NF1_DEgenes.csv", sep = ""), parentId = "syn8506635"), executed = this.file)

