library(data.table)
library(dplyr)
library(synapseClient)

files<-list.files('../data/NF2_only_samples')
files<-gsub("NF2only", "", files)
files<-gsub("mutations.txt", "", files)

panCancerMutationBurden<-function(cancerType){
  mutations<-read.table(file=paste('../data/mutations_by_cancer/',cancerType,'mutations.txt', sep=""), head=TRUE)
  mutn <- mutations[,-1]
  mutations$sums<-rowSums(mutn)
  mutations$ratios<-(mutations$sums/(ncol(mutations)-1))
  mutations$Cancer_Type<-c(rep(cancerType, nrow(mutations)))
  mutations$Hugo_Symbol<-mutations$gene
  mutations<-dplyr::select(mutations, Cancer_Type, Hugo_Symbol, sums, ratios)
  return(mutations)
}

mutationBurden<-lapply(files, function(x){
  print(x)
  df<-panCancerMutationBurden(x)
})

names(mutationBurden)<-files

readMutationFiles<-function(cancerType){
  samples<-read.table(file = paste('../data/NF2_only_samples/NF2only', cancerType, "mutations.txt", sep=""), head = TRUE)
  sample.list<-colnames(samples)[-1]
  sample.list<-head(sample.list, -1)
  sample.list<-gsub("\\.", '-', sample.list)
  
    cancer<-lapply(sample.list, function(i) {
    patientmuts.df <- try(fread(file=paste("../data/tcga_mafs/gdac.broadinstitute.org_",cancerType,".Mutation_Packager_Oncotated_Calls.Level_3.2016012800.0.0/",i,".hg19.oncotator.hugo_entrez_remapped.maf.txt", sep=""), header = TRUE), silent = TRUE)
    patient.id.col <- rep(i, nrow(patientmuts.df))
    try(patientmuts.df$Patient_ID <- patient.id.col)
    return(patientmuts.df)
  })
  cancer.df<-rbindlist(cancer, use.names = TRUE, fill = TRUE)
  try(mutation.count<-dplyr::select(samples, gene, sums))
  try(names(mutation.count)[names(mutation.count) == "sums"] <- "mutation_count_in_NF2")
  try(names(mutation.count)[names(mutation.count) == "gene"] <- "Hugo_Symbol")
  numsamps<-length(sample.list)
  try(mutation.count$mutation_ratio_in_NF2<-(mutation.count$mutation_count_in_NF2/numsamps))
  try(cancer.type.col<-rep(cancerType, nrow(cancer.df)))
  try(cancer.df<- merge(x=cancer.df, y=mutation.count, by = "Hugo_Symbol"))
  try(cancer.df$Cancer_Type <- cancer.type.col)
  try(mutationrates<-mutationBurden[[cancerType]])
  try(cancer.df<-merge(x = cancer.df, y = mutationrates, by = "Hugo_Symbol"))
  return(cancer.df)
}

cancers<- lapply(files, function(x) {
  print(x)
  readMutationFiles(x)
})

pancancer.NF2comuts<-rbindlist(cancers, use.names = TRUE, fill = TRUE)
write.table(pancancer.NF2comuts, file="NF2_Co_Mutated_Genes_in_TCGA.txt", sep = "\t")

pancancer.NF2comuts<-dplyr::select(pancancer.NF2comuts, Cancer_Type.x,Patient_ID,Hugo_Symbol,Chromosome,
                                   Start_position,End_position,Strand,Variant_Classification,Reference_Allele,Tumor_Seq_Allele1,
                                   Tumor_Seq_Allele2,dbSNP_RS,cDNA_Change,Codon_Change,Protein_Change,SwissProt_entry_Id,Description,
                                   sums, mutation_count_in_NF2, ratios)
pancancer.NF2comuts<-dplyr::filter(pancancer.NF2comuts, mutation_count_in_NF2>0)
names(pancancer.NF2comuts)[names(pancancer.NF2comuts) == 'sums'] <- "mutations_in_cancer_type"
names(pancancer.NF2comuts)[names(pancancer.NF2comuts) == 'Cancer_Type.x'] <- "Cancer_Type"
names(pancancer.NF2comuts)[names(pancancer.NF2comuts) == 'ratios'] <- "mut_freq_in_cancer_type"

schema<-synGet("syn8069274")
newcol<-TableColumn(name="mut_freq_in_cancer_type", columnType = "DOUBLE")
schema<-synAddColumn(schema, newcol)
schema<-synStore(schema)

queryResult<-synTableQuery("select * from syn8069274")
queryResult@values["mut_freq_in_cancer_type"] <- pancancer.NF2comuts$mut_freq_in_cancer_type
table<-synStore(queryResult)
