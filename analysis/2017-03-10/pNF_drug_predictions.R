library(synapseClient)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggbeeswarm)

synapseLogin()

mipe.map<-read.table(synGet("syn5522649")@filePath, header = TRUE, sep = ",", quote = "\"") %>%
  select(NCGC.SID, name)
colnames(mipe.map) <- c("Supplier_ID", "Drug")

pNF.screen <- read.table(synGet("syn5637634")@filePath, header = TRUE, sep = "\t", quote = "\"") %>%
  full_join(mipe.map)

pNF.predict <- read.table(synGet("syn8331068")@filePath, header = TRUE, sep = "\t", quote = "\"") %>%
  filter(Hypergeo_pval <= 0.1 & Supplier == "MIPE")

pNF.predict.test <- filter(pNF.screen, pNF.screen$Supplier_ID %in% pNF.predict$Supplier_ID) %>%
  select(Cell, Drug, AUC) %>% 
  spread(Drug, AUC)

rownames(pNF.predict.test) <- pNF.predict.test[,1]
pNF.predict.test <- pNF.predict.test[,-1]

annot <- synTableQuery("SELECT sampleIdentifier, nf1Genotype FROM syn7850572")@values %>%
  distinct() %>%
  slice(1:9) 

names <- annot[,1]
annot <- as.data.frame(annot[,-1])
ip62<-c("+/+")
annot <- rbind(annot, ip62)
rownames(annot) <- c(names, "ipNF06.2A")

##auc data
pheatmap(pNF.predict.test, annotation_row = annot)



##get other response variables from synapse
one <- read.table(synGet("syn5522642")@filePath, header = TRUE, sep = ",", quote = "\"")
unify <- colnames(one)
two <- read.table(synGet("syn5522643")@filePath, header = TRUE, sep = ",", quote = "\"") 
three <- read.table(synGet("syn5522644")@filePath, header = TRUE, sep = ",", quote = "\"") 
four <- read.table(synGet("syn5522645")@filePath, header = TRUE, sep = ",", quote = "\"")
five <- read.table(synGet("syn5522646")@filePath, header = TRUE, sep = ",", quote = "\"")
six <- read.table(synGet("syn5522648")@filePath, header = TRUE, sep = ",", quote = "\"")
seven <-read.table(synGet("syn5522647")@filePath, header = TRUE, sep = ",", quote = "\"") 
eight <- read.table(synGet("syn5522649")@filePath, header = TRUE, sep = ",", quote = "\"")

colnames(two) <- unify
colnames(three) <- unify  
colnames(four) <- unify  
colnames(five) <- unify
colnames(six) <- unify
colnames(seven) <- unify
colnames(eight) <- unify

allDRData <- bind_rows(one, two, three, four, five, six, seven, eight)
names(allDRData$NCGC.SID) <- c("Supplier_ID")

pNF.IC50 <- allDRData %>% 
  filter(allDRData$NCGC.SID %in% pNF.predict$Supplier_ID) %>% 
  select(Cell.line, name, LAC50) %>% 
  spread(name, LAC50)

rownames(pNF.IC50) <- c("ipn02.3", "ipn02.8", "ipNF05.5 (mixed clone)", "ipNF05.5 (single clone)", "ipNF06.2A", "ipNF95.11bC", "ipNF95.6", "ipnNF95.11c")
pNF.IC50 <- pNF.IC50[,-1]

pheatmap(pNF.IC50, annotation_row = annot)


pNF.MR <- allDRData %>% 
  filter(allDRData$NCGC.SID %in% pNF.predict$Supplier_ID) %>% 
  mutate(100-MAXR) %>%
  select(Cell.line, name, `100 - MAXR`) %>%
  spread(name, `100 - MAXR`)

rownames(pNF.MR) <- c("ipn02.3", "ipn02.8", "ipNF05.5 (mixed clone)", "ipNF05.5 (single clone)", "ipNF06.2A", "ipNF95.11bC", "ipNF95.6", "ipnNF95.11c")
pNF.MR <- pNF.MR[,-1]

pheatmap(pNF.MR, annotation_row = annot)


## lm of drug response across cell types 
Cell.line<- c("ipNF02.3 2l", "ipNF02.8", "ipNF05.5 Mixed Clones", "ipNF05.5 Single Clone", 
              "ipnNF95.11C", "ipNF95.11b C/T", "ipNF95.6", "ipNF06.2A")
genotype <- c("+/+", "+/+", "-/-", "-/-", "-/-", "-/-", "-/-", "+/+")

annot<-as.data.frame(cbind(Cell.line, genotype))

filteredData <- allDRData %>% 
  left_join(annot) %>% 
  filter(NCGC.SID %in% pNF.predict$Supplier_ID) %>%
  select(Cell.line, LAC50, MAXR, FAUC, name, genotype) 

lm<-lm(FAUC ~ genotype, filteredData)
summary(lm)

ggplot(filteredData, aes(x = Cell.line, y = MAXR), group = genotype) +
  geom_beeswarm(aes(color = genotype))
  