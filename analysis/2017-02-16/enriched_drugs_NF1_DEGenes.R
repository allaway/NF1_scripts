library(synapseClient)
library(plyr)
library(dplyr)
library(biomaRt)
library(mHG)
library(parallel)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-02-16/enriched_drugs_NF1_DEGenes.R"


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

## make more stringent version of above drug data (increased evidence
## for hitting target, more potency towards target) forget everything
## requiring more than 50uM drug to have an effect, would be very
## challenging to deliver that concentration of drug to target in tumor
drugdat.stringent <- filter(drugdat.filtered, N_quantitative > 1)
drugdat.stringent <-
  filter(drugdat.filtered, MinActivity_nM < 50000)
listofdrugtargets.stri <-
  dlply(.data = drugdat.stringent, .variables = "Structure_ID")

##get NF1 DEgene filenames
files <-
  synQuery("SELECT * from file WHERE parentId=='syn8267685'")$file.id

## core function to test for enriched drug targets (hypergeometric test)
TestForDrugTargets <- function(comut) {
  allcomuts <- unique(comut)
  hyper <- lapply(listofdrugtargets.stri, function(x) {
    N <- length(allcomuts)
    B <- nrow(x$Hugo_Gene)
    lambdas <- as.integer((allcomuts %in% x$Hugo_Gene))
    mHG <- mHG.test(lambdas)$p.value
  })
  
  Structure_ID <- names(listofdrugtargets.stri)
  hypergeo_pval <- t(bind_rows(hyper))
  hyper.df <- as.data.frame(cbind(Structure_ID, hypergeo_pval))
  names(hyper.df) <- c("Structure_ID", "Hypergeo_pval")
  compound.data$Structure_ID <-
    as.character(compound.data$Structure_ID)
  hyper.annot <-
    left_join(hyper.df, compound.data, by = "Structure_ID")
}

cores<-detectCores()
print(cores)

##lapply across all tumor types
hyper <- mclapply(files, function(x) {
  print(x)
  syn <- synGet(x)
  cancer <- syn@fileHandle$fileName
  cancer <- sub("TCGA_", "", cancer)
  cancer <- sub("_NF1_DEgenes.csv", "", cancer)
  
  DEgenes <- read.table(syn@filePath, sep = ",", header = TRUE)
  DEgenes$Gene <- rownames(DEgenes)
  incr.in.NF1 <- filter(DEgenes, logFC < -0.1 & adj.P.Val < 0.1)
  sort.gene <- arrange(DEgenes, logFC)
  
  hyper.annot <- TestForDrugTargets(sort.gene$Gene)
  cancer_type <- rep(cancer, nrow(hyper.annot))
  hyper.annot <- cbind(hyper.annot, cancer_type)
  
}, mc.cores = cores)

##consolidate back into one df, adjust pval for multi corrections, make df with significantly enriched compounds
hyper.df <- bind_rows(hyper)
hyper.df$pval_BHadj <- p.adjust(hyper.df$Hypergeo_pval, method = "BH")
sigs.df <- filter(hyper.df, pval_BHadj < 0.05)
#ids <- count(sigs.df$Structure_ID)

write.table(hyper.df, "NF1_TCGA_DEGenes_enriched_drugs.txt", sep = "\t")
synStore(File("NF1_TCGA_DEGenes_enriched_drugs.txt", parentId = 'syn8292685'), used = files, executed = this.file)
