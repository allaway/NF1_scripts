library(plyr)
library(dplyr)
library(reshape2)
library(synapseClient)
library(ggplot2)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-02-23/PlotSigDEGenesTCGANF1.R"

files <-
  synQuery("SELECT * from file WHERE parentId=='syn8267685'")$file.id

#volcano plots
for(i in files){
  data<-read.table(synGet(i)@filePath, sep = ",", header = TRUE)
  cancer<-synGet(i)@fileHandle$fileName
  cancer<-sub("TCGA_", "", cancer)
  cancer<-sub("_NF1_DEgenes.csv", "", cancer)
  data$gene <- rownames(data)

  ggplot(data = data) +
    geom_point(stat = "identity", aes(x=logFC, y=-log10(adj.P.Val), color = -log10(adj.P.Val)>1.301)) +
    geom_text(aes(x=logFC, y=-log10(adj.P.Val), label = gene), data = data %>% filter(-log10(adj.P.Val)>1.301),
            nudge_x = (0.1), hjust = 0, size = 2, position = ) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    labs(title = paste("NF1 DEgenes in ", cancer, " (NF1 wt vs NF1 mut)", sep = ""), x = "Log10(FC)", y = "-log10(p-val)") +
    theme(plot.title=element_text(hjust=0.5)) +
    guides(color = "none")
  ggsave(paste("sig_NF1_DE_genes_in_",cancer,"_TCGA.png", sep = ""), plot = last_plot(), width = 7, height = 7) 
  synStore(File(paste("sig_NF1_DE_genes_in_",cancer,"_TCGA.png", sep = ""),parentId = "syn8314098"), used = i, executed = this.file) 

}

#waterfall plots
for(i in files){
  data<-read.table(synGet(i)@filePath, sep = ",", header = TRUE)
  cancer<-synGet(i)@fileHandle$fileName
  cancer<-sub("TCGA_", "", cancer)
  cancer<-sub("_NF1_DEgenes.csv", "", cancer)
  data$gene <- rownames(data)
  data$gene<-factor(data$gene, levels = data$gene[order(-data$logFC)])
  if(nrow(filter(data, -log10(adj.P.Val)>1.301))>0) {
  ggplot(data = data %>% filter(-log10(adj.P.Val)>1.301)) +
    geom_bar(stat = "identity", aes(x=gene, y=-logFC, fill = -logFC)) +    
    scale_fill_gradient2(low="darkblue", high="red", midpoint=0) +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
    labs(title = paste("NF1 DEgenes in ", cancer, " (NF1 mut vs NF1 wt), adj. p<0.1", sep = ""), y = "Log10(FC)", x = "Gene") +
    theme(plot.title=element_text(hjust=0.5)) +
    guides(color = "none") 
  ggsave(paste("sig_NF1_DE_genes_in_",cancer,"_TCGA_waterfall.png", sep = ""), plot = last_plot(), width = 7, height = 7) 
  synStore(File(paste("sig_NF1_DE_genes_in_",cancer,"_TCGA_waterfall.png", sep = ""),parentId = "syn8384412"), used = i, executed = this.file) 
  } else {
    print("no sig DEgenes")
  }
}
