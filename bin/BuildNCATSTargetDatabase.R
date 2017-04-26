library(synapseClient)
library(plyr)
library(dplyr)
library(biomaRt)
synapseLogin()

this.file = ""


## pull drug data and filter for human targets, and eliminate drugs with
## 0 quantitative effects measured
drugdat <- synTableQuery("SELECT * FROM syn7341038")
drugdat <- as.data.frame(drugdat@values)
drugdat.filtered <- filter(drugdat, Organism == "Homo sapiens")
drugdat.filtered <- filter(drugdat.filtered, N_quantitative != 0)

## obtain data to map Uniprot id with Hugo Genes
mart <- useMart("ENSEMBL_MART_ENSEMBL", "hsapiens_gene_ensembl")
bm <-
  getBM(attributes = c("hgnc_symbol", "uniprot_swissprot"),
        mart = mart)
bm <- filter(bm, uniprot_swissprot != "")
colnames(bm) <- c("Hugo_Gene", "Uniprot_accession_numbers")

## map uniprot targets to hugo genes - generate list of dataframes, each
## data frame lists targets of drug
drugdat.filtered$Uniprot_accession_numbers <-
  sub(",.", "", drugdat.filtered$Uniprot_accession_numbers)
drugdat.filtered <-
  left_join(drugdat.filtered, bm, by = "Uniprot_accession_numbers")
listofdrugtargets <-
  dlply(.data = drugdat.filtered, .variables = "Structure_ID")

## some uniprot ids do not successfully map to hugo genes, note them
## here, these will not show up in analysis
unannotated <-
  filter(drugdat.filtered, is.na(drugdat.filtered$Hugo_Gene))

## pull evotec data to map compound names to structure IDs
compound.data <- synTableQuery("SELECT * FROM syn8118065")
compound.data <- as.data.frame(compound.data@values)
compound.data <-
  compound.data[!duplicated(compound.data$Structure_ID),]

## screening data NF1/NF2
NF2.ncats <- read.table(synGet("syn8314523")@filePath, header = TRUE, comment.char = "")
NF1.ncats
