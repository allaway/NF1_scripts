library(synapseClient)
library(ggplot2)
library(dplyr)
library(rJava)
library(rcdk)
synapseLogin()
library(pheatmap)

x<-synGet("syn8295452")@filePath
data<-read.table(x)
sigs<-filter(data, pval_BHadj<=0.1)
sig.data<-select(sigs, Structure_ID, Original_molecule_SMILES, Supplier_Data_1,Supplier_Data_2,Supplier_Data_3)
sig.data<-distinct(sig.data)

mol<-parse.smiles(as.character(sig.data$Original_molecule_SMILES))

#write.molecules(mol, filename = 'NF1_TCGA_drugs.sdf')
#view.molecule.2d(mol)

fps <- lapply(mol, get.fingerprint, type="extended")
fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
rownames(fp.sim) <- sig.data$Structure_ID
colnames(fp.sim) <- sig.data$Structure_ID
fp.dist <- 1- fp.sim
pheatmap(fp.dist)

clust<-hclust(as.dist(fp.dist))
plot(clust)

