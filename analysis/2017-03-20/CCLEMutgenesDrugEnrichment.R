library(synapseClient)
library(reshape2)
library(plyr)
library(dplyr)
library(biomaRt)
library(mHG)
library(parallel)
synapseLogin()

this.file = "https://raw.githubusercontent.com/allaway/NF1_scripts/master/analysis/2017-03-20/CCLEMutgenesDrugEnrichment.R"

mutations<-read.table(synGet("syn7466552")@filePath)
gene<-rownames(mutations)
mutations<-cbind(gene, mutations)

nf1muts<-as.data.frame(t(filter(mutations, gene == "NF1")))
nf1muts$samps<-rownames(nf1muts)
nf1muts<-filter(nf1muts, V1 == "1")$samps
nf1mut<-dplyr::select(mutations, one_of(nf1muts))
nf1mut$sums <- rowSums(nf1mut)
nf1mut$gene <- rownames(nf1mut)
nf1mut <- nf1mut %>% 
  dplyr::select(gene, sums) %>%
  arrange(desc(sums))

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

print(detectCores())

## core function to test for enriched drug targets (hypergeometric test)
TestForDrugTargets <- function(comut) {
  allcomuts <- unique(comut)
  hyper <- mclapply(listofdrugtargets.stri, function(x) {
    N <- length(allcomuts)
    B <- nrow(x$Hugo_Gene)
    lambdas <- as.integer((allcomuts %in% x$Hugo_Gene))
    mHG <- mHG.test(lambdas)$p.value
  }, mc.cores=detectCores())
  
  Structure_ID <- names(listofdrugtargets.stri)
  hypergeo_pval <- t(bind_rows(hyper))
  hyper.df <- as.data.frame(cbind(Structure_ID, hypergeo_pval))
  names(hyper.df) <- c("Structure_ID", "Hypergeo_pval")
  compound.data$Structure_ID <-
    as.character(compound.data$Structure_ID)
  hyper.annot <-
    left_join(hyper.df, compound.data, by = "Structure_ID")
}

testgenes<-unique(nf1mut$gene)

hyper<-TestForDrugTargets(testgenes)

write.table(hyper, "CCLE_mutation_enriched_drugs.txt", sep = "\t")
synStore(File("CCLE_mutation_enriched_drugs.txt", parentId = "syn8496331"), used = c("syn7466552","syn7341038", "syn8118065"), executed = this.file)


##try again with mHG NF1 comuts - different list - ranked by hypergeo pva
mutations<-read.table(synGet("syn8495196")@filePath) %>%
  arrange(p_value)
testgenes<-unique(mutations$gene)
hyper<-TestForDrugTargets(testgenes)

write.table(hyper, "CCLE_SIGmutation_enriched_drugs.txt", sep = "\t")
synStore(File("CCLE_SIGmutation_enriched_drugs.txt", parentId = "syn8496331"), used = c("syn8495196","syn7341038", "syn8118065"), executed = this.file)


