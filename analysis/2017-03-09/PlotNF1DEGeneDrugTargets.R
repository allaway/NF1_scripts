library(synapseClient)
library(biomaRt)
library(dplyr)
library(Matrix)
library(pheatmap)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-03-09/PlotNF1DEGeneDrugTargets.R"

de.drugs<-read.table(synGet("syn8295452")@filePath, na.strings=c("", "NA")) 
de.drugs <- de.drugs %>% filter(Hypergeo_pval<=0.01)

##get targets using EVOTEC IDs and biomart
targets <- synTableQuery("SELECT * FROM syn7341038")
targets <- as.data.frame(targets@values)
targets <- filter(targets, Organism == 'Homo sapiens')

## obtain data to map Uniprot id with Hugo Genes
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
bm <- getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"),
            mart = mart)
bm <- filter(bm, uniprot_swissprot != "")
colnames(bm) <- c("Hugo_Gene", "Uniprot_accession_numbers")

## map uniprot targets to hugo genes - generate list of dataframes, each  
## data frame lists targets of drug
targets <- left_join(targets, bm, by = "Uniprot_accession_numbers")
targets <- dplyr::select(targets, Structure_ID, Protein_names, Hugo_Gene)

targets<-left_join(de.drugs, targets)
cancers<-unique(targets$cancer_type)

targets$Hugo_Gene<-as.factor(targets$Hugo_Gene)


##plot top 50 targets in each cancer
for(i in cancers){
  foo <- filter(targets, cancer_type == i & !is.na(Hugo_Gene))
  bar <- as.data.frame(summary(foo$Hugo_Gene, maxsum = 20000)[-1])
  names(bar) <- "count"
  bar$count <- as.numeric(bar$count)
  bar$gene <- rownames(bar)
  bar <- bar %>% filter(!is.na(gene)) %>%
    arrange(desc(count)) %>%
    slice(1:50)
    
  ggplot(bar, aes(x=gene, y=count)) + 
    geom_bar(stat="identity") +
    scale_x_discrete(limits = bar$gene) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  #ggsave(paste("top_50_targets_in_", i, ".png", sep = ""))
  #synStore(File(paste("top_50_targets_in_", i, ".png", sep = ""), parentId = "syn8404515"), used = c("syn8295452", "syn7341038"), executed=this.file)
}


##generate heatmap to cluster molecules for cancers
for(i in cancers){
  foo <- filter(targets, cancer_type == i & !is.na(Hugo_Gene) & Hugo_Gene != "")
  foo <- droplevels(foo)
  idx<-cbind(foo$Structure_ID, as.numeric(as.factor(foo$Structure_ID)), as.data.frame(foo$Hugo_Gene), as.numeric(as.factor(foo$Hugo_Gene)))
  colnames(idx) <- c("Structure_ID", "i", "Hugo_Gene", "j")
  mat <- 1*as.matrix(sparseMatrix(i = idx$i, j = idx$j, dimnames = list(unique(idx$Structure_ID), (unique(idx$Hugo_Gene[order(idx$Hugo_Gene)])))))
  mat2<-as.matrix(sparseMatrix(i = idx$i, j = idx$j, dimnames = list(unique(idx$Structure_ID), (unique(idx$Hugo_Gene[order(idx$Hugo_Gene)])))))
  image(mat2)
  #png(paste("compound_relationship_heatmap_DEgenes_,",i,".png", sep = ""))
  pheatmap(mat)
  #dev.off()
  #synStore(File(paste("compound_relationship_heatmap_DEgenes_,",i,".png", sep = ""), parentId = "syn8404515"), used = c("syn8295452", "syn7341038"), executed=this.file)
}


##generate heatmap to cluster molecules for cancers, adjusted by logFC of each gene
files <-
  synQuery("SELECT * from file WHERE parentId=='syn8267685'")$file.id
x <- "syn8269290"
for(x in files){
  print(x)
  syn <- synGet(x)
  i <- syn@fileHandle$fileName
  i <- sub("TCGA_", "", i)
  i <- sub("_NF1_DEgenes.csv", "", i)
    
  DEgenes <- read.table(syn@filePath, sep = ",", header = TRUE) 
  DEgenes$gene <- rownames(DEgenes) 
  
  foo <- filter(targets, cancer_type == i & !is.na(Hugo_Gene) & Hugo_Gene != "")
  foo <- droplevels(foo)
  idx<-cbind(foo$Structure_ID, as.numeric(as.factor(foo$Structure_ID)), as.data.frame(foo$Hugo_Gene), as.numeric(as.factor(foo$Hugo_Gene)))
  colnames(idx) <- c("Structure_ID", "i", "Hugo_Gene", "j")
  mat <- as.data.frame(1*as.matrix(sparseMatrix(i = idx$i, j = idx$j, dimnames = list(unique(idx$Structure_ID), (unique(idx$Hugo_Gene[order(idx$Hugo_Gene)]))))))
  mat <- dplyr::select(mat, one_of(DEgenes$gene))
  DEgenes <- dplyr::select(DEgenes, logFC)
  for(k in colnames(mat)) {
    mat[,k]<-mat[,k]*DEgenes[k,1]
  }
  mat<-as.matrix(mat)
  #mat2<-as.matrix(sparseMatrix(i = idx$i, j = idx$j, dimnames = list(unique(idx$Structure_ID), (unique(idx$Hugo_Gene[order(idx$Hugo_Gene)])))))
  #image(mat2)
  #png(paste("LogFC_adjusted_compound_relationship_heatmap_DEgenes_,",i,".png", sep = ""))
  res<-pheatmap(mat)
  res
  #dev.off()
  #synStore(File(paste("LogFC_adjusted_compound_relationship_heatmap_DEgenes_,",i,".png", sep = ""), parentId = "syn8404515"), used = c("syn8295452", "syn7341038", syn), executed=this.file)
  
  clust.r<- cutree(res$tree_row, k=10)
  clust.c<- cutree(res$tree_col, k=10)
  clust.r
  clust.c
}

 

