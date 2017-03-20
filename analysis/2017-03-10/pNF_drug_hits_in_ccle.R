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

mipe.map<-read.table(synGet("syn5522649")@filePath, header = TRUE, sep = ",", quote = "\"") %>%
  select(NCGC.SID, name)
colnames(mipe.map) <- c("Supplier_ID", "Drug")

ccle.predictions<-read.table(synGet("syn8496371")@filePath) %>%
  filter(Hypergeo_pval <= 0.05 & Supplier == "MIPE") %>%
  left_join(mipe.map)

##manually plot example drugs
ggplot(data = ccle.drug %>% select(NVP.TAE684, line, NF1), aes(x=NF1, y=NVP.TAE684)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(lapatinib, line, NF1), aes(x=NF1, y=lapatinib)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(saracatinib, line, NF1), aes(x=NF1, y=saracatinib)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(GSK.3.inhibitor.IX, line, NF1), aes(x=NF1, y=GSK.3.inhibitor.IX)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(linsitinib, line, NF1), aes(x=NF1, y=linsitinib)) + 
  geom_boxplot(aes(fill = NF1))

##tested vs sig comuts now
ccle.predictions2<-read.table(synGet("syn8497599")@filePath) %>%
  filter(Hypergeo_pval <= 0.05 & Supplier == "MIPE") %>%
  left_join(mipe.map)

##manually plot example drugs
ggplot(data = ccle.drug %>% select(NVP.TAE684, line, NF1), aes(x=NF1, y=NVP.TAE684)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(lapatinib, line, NF1), aes(x=NF1, y=lapatinib)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(saracatinib, line, NF1), aes(x=NF1, y=saracatinib)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(GSK.3.inhibitor.IX, line, NF1), aes(x=NF1, y=GSK.3.inhibitor.IX)) + 
  geom_boxplot(aes(fill = NF1))

ggplot(data = ccle.drug %>% select(linsitinib, line, NF1), aes(x=NF1, y=linsitinib)) + 
  geom_boxplot(aes(fill = NF1))
