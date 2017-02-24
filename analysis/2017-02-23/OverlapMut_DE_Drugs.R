library(synapseClient)
library(dplyr)
library(rJava)
library(rcdk)
synapseLogin()
library(pheatmap)
library(viridis)
library(ape)

DE<-synGet("syn8295452")@filePath
DE.df<-read.table(DE)
DE.df<-filter(DE.df, Hypergeo_pval<=0.005)
dedrugs<-unique(DE.df$Structure_ID)

coMut<-synGet("syn8292948")@filePath
coMut.df<-read.table(coMut)
coMut.df<-filter(coMut.df, Hypergeo_pval<=0.005)
comutdrugs<-unique(coMut.df$Structure_ID)

overlap<-intersect(dedrugs, comutdrugs)

smiledat<-unique(dplyr::select(DE.df, Structure_ID, Original_molecule_SMILES))
smiledat<-filter(smiledat, Structure_ID %in% overlap)

mol<-parse.smiles(as.character(smiledat$Original_molecule_SMILES))

#write.molecules(mol, filename = 'NF1_TCGA_drugs.sdf')
#view.molecule.2d(mol

fps <- lapply(mol, get.fingerprint, type="extended")
fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
rownames(fp.sim) <- smiledat$Structure_ID
colnames(fp.sim) <- smiledat$Structure_ID
fp.dist <- 1- fp.sim
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000))

clust<-hclust(as.dist(fp.dist))
phylo<-as.phylo(clust)
plot(phylo, type = "fan", cex = 0.65, label.offset = 0.01)
