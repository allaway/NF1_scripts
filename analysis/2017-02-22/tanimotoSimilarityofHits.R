library(synapseClient)
library(ggplot2)
library(dplyr)
library(rJava)
library(rcdk)
synapseLogin()
library(pheatmap)
library(viridis)
library(ape)

this.file = "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-02-22/tanimotoSimilarityofHits.R"

###NF1 TCGA DEGene Drug Hits 
x<-synGet("syn8295452")@filePath
data<-read.table(x)
sigs<-filter(data, Hypergeo_pval<=0.0005)
sig.data<-dplyr::select(sigs, Structure_ID, Original_molecule_SMILES, Supplier_Data_1,Supplier_Data_2,Supplier_Data_3)
sig.data<-distinct(sig.data)

mol<-parse.smiles(as.character(sig.data$Original_molecule_SMILES))

#write.molecules(mol, filename = 'NF1_TCGA_drugs.sdf')
#view.molecule.2d(mol)

fps <- lapply(mol, get.fingerprint, type="extended")
fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
rownames(fp.sim) <- sig.data$Structure_ID
colnames(fp.sim) <- sig.data$Structure_ID
fp.dist <- 1- fp.sim

svg("NF1_DEgenes_enriched_drugs_tanimoto_heatmap.svg")
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000) , show_colnames = FALSE, show_rownames = FALSE)
dev.off()

synStore(File("NF1_DEgenes_enriched_drugs_tanimoto_heatmap.svg", parentId = "syn8327212"), used = "syn8295452", executed = this.file)

png("NF1_DEgenes_enriched_drugs_tanimoto_heatmap.png")
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000), show_colnames = FALSE, show_rownames = FALSE)
dev.off()

synStore(File("NF1_DEgenes_enriched_drugs_tanimoto_heatmap.png", parentId = "syn8327212"), used = "syn8295452", executed = this.file)

clust<-hclust(as.dist(fp.dist))
png("NF1_DEgenes_enriched_drugs_dendro.png")
plot(clust)
dev.off()

synStore(File("NF1_DEgenes_enriched_drugs_dendro.png", parentId = "syn8327212"), used = "syn8295452", executed = this.file)

svg("NF1_DEgenes_enriched_drugs_dendro.svg")
plot(clust)
dev.off()

synStore(File("NF1_DEgenes_enriched_drugs_dendro.svg", parentId = "syn8327212"), used = "syn8295452", executed = this.file)

phylo<-as.phylo(clust)
png("NF1_DEgenes_enriched_drugs_dendro_circular.png")
plot.phylo(phylo, type = "fan", cex = 0.65, label.offset = 0.01, tip.color = "black", edge.color = rainbow(n = 2*nrow(fp.dist)))
dev.off()

synStore(File("NF1_DEgenes_enriched_drugs_dendro_circular.png", parentId = "syn8327212"), used = "syn8295452", executed = this.file)


svg("NF1_DEgenes_enriched_drugs_dendro_circular.svg")
plot.phylo(phylo, type = "fan", cex = 0.65, label.offset = 0.01, tip.color = "black", edge.color = rainbow(n = 2*nrow(fp.dist)))
dev.off()

synStore(File("NF1_DEgenes_enriched_drugs_dendro_circular.svg", parentId = "syn8327212"), used = "syn8295452", executed = this.file)


fit<-as.data.frame(cmdscale(fp.dist))
fit$structure_id<-rownames(fit)
colnames(fit) <- c("dimension_1", "dimension_2", "structure_id")
ggplot(data = fit, aes(x=dimension_1, y=dimension_2)) +
  geom_point(stat = "identity", color = "blue") +
  geom_text(aes(label=structure_id), vjust = -1.5, size = 2) +
  coord_fixed(ratio = 1)

ggsave("MDS_of_NF1_DE_enriched_drugs.png", height = 5, width = 5)

synStore(File("MDS_of_NF1_DE_enriched_drugs.png", parentId = "syn8327212"), used = "syn8295452", executed = this.file)


###NF1 TCGA Comutant Hits
x<-synGet("syn8292948")@filePath
data<-read.table(x)
sigs<-filter(data, Hypergeo_pval<=0.0005)
sig.data<-dplyr::select(sigs, Structure_ID, Original_molecule_SMILES, Supplier_Data_1,Supplier_Data_2,Supplier_Data_3)
sig.data<-distinct(sig.data)

mol<-parse.smiles(as.character(sig.data$Original_molecule_SMILES))

#write.molecules(mol, filename = 'NF1_TCGA_drugs.sdf')
#view.molecule.2d(mol)

fps <- lapply(mol, get.fingerprint, type="extended")
fp.sim <- fp.sim.matrix(fps, method = "tanimoto")
rownames(fp.sim) <- sig.data$Structure_ID
colnames(fp.sim) <- sig.data$Structure_ID
fp.dist <- 1- fp.sim

svg("NF1_COMUTgenes_enriched_drugs_tanimoto_heatmap.svg")
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000), main = "Tanimoto Similarity", show_colnames = FALSE, show_rownames = FALSE)
dev.off()

synStore(File("NF1_COMUTgenes_enriched_drugs_tanimoto_heatmap.svg", parentId = "syn8327214"), used = "syn8292948", executed = this.file)

png("NF1_COMUTgenes_enriched_drugs_tanimoto_heatmap.png")
pheatmap(fp.sim, border_color = NA, color = magma(n = 10000), main = "Tanimoto Similarity", show_colnames = FALSE, show_rownames = FALSE)
dev.off()

synStore(File("NF1_COMUTgenes_enriched_drugs_tanimoto_heatmap.png", parentId = "syn8327214"), used = "syn8292948", executed = this.file)


clust<-hclust(as.dist(fp.dist))

png("NF1_COMUTgenes_enriched_drugs_dendro.png")
plot(clust)
dev.off()

synStore(File("NF1_COMUTgenes_enriched_drugs_dendro.png", parentId = "syn8327214"), used = "syn8292948", executed = this.file)

svg("NF1_COMUTgenes_enriched_drugs_dendro.svg")
plot(clust)
dev.off()

synStore(File("NF1_COMUTgenes_enriched_drugs_dendro.svg", parentId = "syn8327214"), used = "syn8292948", executed = this.file)


phylo<-as.phylo(clust)

png("NF1_COMUTgenes_enriched_drugs_dendro_circular.png")
plot.phylo(phylo, type = "fan", cex = 0.65, label.offset = 0.01, tip.color = "black", edge.color = rainbow(n = 2*nrow(fp.dist)))
dev.off()

synStore(File("NF1_COMUTgenes_enriched_drugs_dendro_circular.png", parentId = "syn8327214"), used = "syn8292948", executed = this.file)

svg("NF1_COMUTgenes_enriched_drugs_dendro_circular.svg")
plot.phylo(phylo, type = "fan", cex = 0.65, label.offset = 0.01, tip.color = "black", edge.color = rainbow(n = 2*nrow(fp.dist)))
dev.off()

synStore(File("NF1_COMUTgenes_enriched_drugs_dendro_circular.svg", parentId = "syn8327214"), used = "syn8292948", executed = this.file)


fit<-as.data.frame(cmdscale(fp.dist))
fit$structure_id<-rownames(fit)
colnames(fit) <- c("dimension_1", "dimension_2", "structure_id")
ggplot(data = fit, aes(x=dimension_1, y=dimension_2)) +
  geom_point(stat = "identity", color = "blue") +
  geom_text(aes(label=structure_id), vjust = -1.5, size = 2) +
  coord_fixed()
ggsave("MDS_of_NF1_CoMut_enriched_drugs.png", height = 5, width = 5)

synStore(File("MDS_of_NF1_CoMut_enriched_drugs.png", parentId = "syn8327214"), used = "syn8292948", executed = this.file)

