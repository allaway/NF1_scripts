library(synapseClient)
library(dplyr)
library(ggplot2)
synapseLogin()

ccle.mut<-as.data.frame(read.table(synGet("syn7466552")@filePath))
ccle.mut[,"gene"] <- rownames(ccle.mut)
NFmut <- as.data.frame(t(filter(ccle.mut, gene == "NF1")))
NFmut[,"line"]<- rownames(NFmut)
colnames(NFmut) <- c("NF1", "line")
#NFmut <- NFmut[, NFmut==1
#mutcells <- colnames(NFmut)

ccle.drug<-read.table(synGet("syn7466611")@filePath)
ccle.drug[,"line"] <- rownames(ccle.drug)
ccle.drug <- ccle.drug %>% 
  filter(line %in% colnames(ccle.mut)) %>% 
  left_join(NFmut)

ccle.predictions<-read.table(synGet("syn8496371")@filePath) %>%
  filter(Hypergeo_pval <= 0.05 & 


ggplot(data = ccle.drug %>% select(NVP.TAE684, line, NF1), aes(x=NF1, y=NVP.TAE684)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(, line, NF1), aes(x=NF1, y=)) + 
  geom_boxplot(aes(fill = NF1))
