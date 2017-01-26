library(synapseClient)
library(data.table)
library(dplyr)

##get project genie mutation data
geniemuts<-synGet("syn7851250")@filePath
muts<-read.table(file=geniemuts, header = TRUE, sep = "\t")

##get clinical annotations
clinannot<-synGet('syn7851246')@filePath
clinannots<-read.table(file=clinannot, header = TRUE, sep = "\t")
cancertype<-select(clinannots, SAMPLE_ID, ONCOTREE_CODE)
colnames(cancertype)[colnames(cancertype)=='SAMPLE_ID'] <- "Tumor_Sample_Barcode"

muts<-merge(muts, cancertype, by="Tumor_Sample_Barcode")
muts<-select(muts, ONCOTREE_CODE, Tumor_Sample_Barcode, Hugo_Symbol, Chromosome, Start_Position, End_Position, Strand, 
             Variant_Classification, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, dbSNP_RS, HGVSc, HGVSp)

##determine number of mutations for all genes by cancer type
getMutCount<-function(mutationframe, cancertype) {
  foo<-filter(mutationframe, ONCOTREE_CODE==cancertype)
  foo<-select(foo, ONCOTREE_CODE, Tumor_Sample_Barcode, Hugo_Symbol)
  foo<-distinct(foo)
  df<-count(foo, Hugo_Symbol)
  df$mut_freq_in_cancer_type<-(df$n/length(unique(foo$Tumor_Sample_Barcode)))
  return(df)
}
mutationrates<-lapply(unique(muts$ONCOTREE_CODE), function(x) {
  print(x)
  getMutCount(muts, x)
})

names(mutationrates)<-unique(muts$ONCOTREE_CODE)

for(i in names(mutationrates)){
  mutationrates[[i]]$ONCOTREE_CODE<-rep(i, nrow(mutationrates[[i]]))
}

mutationrates<-rbindlist(mutationrates)

##find NF1 and NF2 mutations in dataset
NF1<-filter(muts, Hugo_Symbol == "NF1")
NF2<-filter(muts, Hugo_Symbol == "NF2")

##get tumor sample barcodes from NF1 mutant tumors
NF1samps<-NF1$Tumor_Sample_Barcode
NF2samps<-NF2$Tumor_Sample_Barcode

##filter all project genie data by these barcodes
NF1mut.samps<-filter(muts, Tumor_Sample_Barcode %in% NF1samps)
NF2mut.samps<-filter(muts, Tumor_Sample_Barcode %in% NF2samps)

##determine mutation rate for all 

##select gene and tumor ID, and reduce observations to account for multiple mutations within one gene in a sample
##to determine mutation frequency
NF1mut.reduced<-distinct(select(NF1mut.samps, Hugo_Symbol, Tumor_Sample_Barcode, ONCOTREE_CODE))
NF2mut.reduced<-distinct(select(NF2mut.samps, Hugo_Symbol, Tumor_Sample_Barcode, ONCOTREE_CODE))

mutationrates.inNF1<-lapply(unique(NF1mut.reduced$ONCOTREE_CODE), function(x) {
  print(x)
  getMutCount(NF1mut.reduced, x)
})

names(mutationrates.inNF1)<-unique(NF1mut.reduced$ONCOTREE_CODE)

for(i in names(mutationrates.inNF1)){
  mutationrates.inNF1[[i]]$ONCOTREE_CODE<-rep(i, nrow(mutationrates.inNF1[[i]]))
}

mutationrates.inNF1<-rbindlist(mutationrates.inNF1)
colnames(mutationrates.inNF1)[colnames(mutationrates.inNF1)=='n'] <- "mut_count_in_nf1_mut"
mutationrates.inNF1<-mutationrates.inNF1[,-3]

mutationrates.inNF2<-lapply(unique(NF2mut.reduced$ONCOTREE_CODE), function(x) {
  print(x)
  getMutCount(NF2mut.reduced, x)
})

names(mutationrates.inNF2)<-unique(NF2mut.reduced$ONCOTREE_CODE)

for(i in names(mutationrates.inNF2)){
  mutationrates.inNF2[[i]]$ONCOTREE_CODE<-rep(i, nrow(mutationrates.inNF2[[i]]))
}

mutationrates.inNF2<-rbindlist(mutationrates.inNF2)
colnames(mutationrates.inNF2)[colnames(mutationrates.inNF2)=='n'] <- "mut_count_in_NF2_mut"
mutationrates.inNF2<-mutationrates.inNF2[,-3]

NF1mut.reduced<-merge(NF1mut.samps,mutationrates.inNF1, by=c("Hugo_Symbol", "ONCOTREE_CODE"))
NF2mut.reduced<-merge(NF2mut.samps,mutationrates.inNF1, by=c("Hugo_Symbol", "ONCOTREE_CODE"))

NF1mut.final<-merge(x=NF1mut.reduced, y=mutationrates, by=c("Hugo_Symbol", "ONCOTREE_CODE"))
colnames(NF1mut.final)[colnames(NF1mut.final)=='n'] <- "mut_count_in_cancer_type"
NF1mut.final.trim<-apply(NF1mut.final, 1:2, function(x) strtrim(x, 1000))
NF1mut.final.trim<-as.data.frame(NF1mut.final.trim)

NF2mut.final<-merge(x=NF2mut.reduced, y=mutationrates, by=c("Hugo_Symbol", "ONCOTREE_CODE"))
colnames(NF2mut.final)[colnames(NF2mut.final)=='n'] <- "mut_count_in_cancer_type"

makeTable <- function(df,tableName,projectId){
  library(synapseClient)
  synapseLogin()
  tcresult<-as.tableColumns(df)
  cols<-tcresult$tableColumns
  getType = sapply(cols,function(x) return(x@columnType))
  fileHandleId<-tcresult$fileHandleId
  project<- synGet(projectId)
  schema<-TableSchema(name=tableName, parent=project, columns=cols)
  table<-Table(schema, fileHandleId)
  table<-synStore(table, retrieveData=TRUE)
}

makeTable(NF1mut.final.trim, "NF1 Co-mutated Genes in Genie", "syn7154892")
makeTable(NF2mut.final, "NF2 Co-mutated Genes in Genie", "syn7154892")


