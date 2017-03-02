library(synapseClient)
library(dplyr)
library(webchem)
library(ChemmineR)
library(ChemmineOB)
synapseLogin()

x<-synTableQuery("SELECT * FROM syn8118065")@values

mols<-unique(x$Original_molecule_SMILES)

token = '4dcc6a4b-7e9a-4ba8-9678-3f6a5207ca5d'

##map as many smiles as possible to CSID
csid<-lapply(mols, function(j){
  y<-get_csid(j, token = token)
  return(c(j,y))
})

csid2<-as.data.frame(csid)
csid2$SMILES <- mols
rownames(csid2) <- c()

unmapped<-filter(csid2, is.na(csid))
mols_unm<-unique(unmapped$SMILES)

##use to map to CID: https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
write.table(as.data.frame(mols), "allmols.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
cids<-read.table(synGet("syn8361218")@filePath, fill = NA, comment.char = "")

colnames(cids) <- c("SMILES", "CID")
allids<-full_join(cids, csid2)
names(allids)<-c("Original_molecule_SMILES", "CID", "CSID")
test<-full_join(allids,x)

write.table(test, "cid_csid_map.txt", sep = "\t")

_##map smiles to 'names' using above resource
write.table(as.data.frame(test$CID), "forsynonyms.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
synonyms<-read.table(synGet("syn8361358")@filePath, quote=NULL, comment='', sep = "\t")

cidcts<-sapply(allids$CSID, function(k){
  pic_prop(k, properties = "names", first = TRUE)
})
 b
