##create a script to parse cBIOportal data from github @sgosline/RASPathwaySig
#perhaps this will facilitate analysis across the different datasets...
library(cgdsr)
library(data.table)

script.dir <- dirname(sys.frame(1)$ofile)
all.genes<-unique(fread('../../data/ucsc_kgXref_hg19_2015_10_29.csv')$geneSymbol)

#'getSamplesForDisease creates a unified mapping of all samples
#'to various cell lines and disease profiles so that when
#'all are joined we can compare one to another
getSamplesForDisease<-function(dis='',study='tcga'){
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  all.studies<-getCancerStudies(mycgds)
  ind=grep(paste(tolower(dis),paste(study,'$',sep=''),sep='_'),all.studies$cancer_study_id)
  if(length(ind)==0)
    return(c())
  mycancerstudy<-all.studies$cancer_study_id[ind]
  
  sampList<-lapply(mycancerstudy,function(cs){
    caseLists<-getCaseLists(mycgds,cs)
    samps<-unlist(strsplit(caseLists[match(paste(cs,'all',sep='_'),caseLists[,1]),5],split=' '))
    #      allprofs<-getGeneticProfiles(mycgds,cs)[,1]
    print(paste('Found',length(samps),'for',cs))
    return(samps)
  })
  all.samps<-unique(unlist(sampList))
  print(paste("Found",length(all.samps),'samples for',study,dis))
  return(all.samps)
}

#various disease types in cbioporta.
broad.cancer.types=c('brca','cellline','lcll','desm','dlbc','esca','hnsc','luad','mbl','skcm','mm','nsclc','es','prad')
mskcc.cancer.types=c('acyc','acbc','blca','coadread','luad','mpnst','thyroid','prad','hnc','sarc','scco')
tcga.cancer.types<-c('laml','acc','blca','lgg','brca','cesc','chol','coadread','esca','gbm','hnsc','kich','kirc','kirp','lihc','luad','lusc','dlbc','lgggbm','ov','nsclc','paad','thca','pcpg','prad','sarc','skcm','stad','tgct','thym','ucs','ucec','uvm')#meso has no sequence data

##not all have counts
cell.line.tiss<-c('CENTRAL_NERVOUS_SYSTEM','BONE','PROSTATE','STOMACH','URINARY_TRACT','OVARY','HAEMATOPOIETIC_AND_LYMPHOID_TISSUE','KIDNEY','THYROID','SKIN','SOFT_TISSUE','SALIVARY_GLAND','LUNG','PLEURA','LIVER','ENDOMETRIUM','PANCREAS','BREAST','UPPER_AERODIGESTIVE_TRACT','LARGE_INTESTINE','AUTONOMIC_GANGLIA','OESOPHAGUS','BILIARY_TRACT','SMALL_INTESTINE')

getDisMutationData<-function(dis='',study='tcga'){
  
  mycgds = CGDS("http://www.cbioportal.org/public-portal/")
  all.studies<-getCancerStudies(mycgds)
  
  if(tolower(dis)=='alltcga')
    dval=''
  else
    dval=dis
  #if disease is blank will get all diseases
  ind=grep(paste(tolower(dval),paste(study,'$',sep=''),sep='_'),all.studies$cancer_study_id)
  print(paste('found',length(ind),study,'samples for disease',dis))
  
  if(length(ind)==0)
    return(NULL)
  
  mycancerstudy<-all.studies$cancer_study_id[ind]
  
  if(study=='tcga')
    mycancerstudy=intersect(mycancerstudy,paste(tcga.cancer.types,study,sep='_'))
  else if(study=='mskcc')
    mycancerstudy=intersect(mycancerstudy,paste(mskcc.cancer.types,study,sep='_'))
  else if(study=='broad')
    mycancerstudy=intersect(mycancerstudy,paste(broad.cancer.types,study,sep='_'))
  
  expr.list<-lapply(mycancerstudy,function(cs){
    print(paste(cs,study,'Mutation data'))
    caseLists<-getCaseLists(mycgds,cs)
    allprofs<-getGeneticProfiles(mycgds,cs)[,1]
    profile=allprofs[grep('mutations',allprofs)]
    seqSamps=caseLists$case_list_id[grep('sequenced',caseLists$case_list_id)]
    gene.groups=split(all.genes, ceiling(seq_along(all.genes)/400))
    dat<-lapply(gene.groups,function(g) getProfileData(mycgds,g,profile,seqSamps))
    ddat<-matrix()
    for(i in which(sapply(dat,nrow)!=0)){
      ddat<-cbind(ddat,dat[[i]])
    }
    nans<-which(apply(ddat,2,function(x) all(is.nan(x)||is.na(x))))
    # nas<-which(apply(ddat,2,function(x) all(is.na(x))))
    ddat<-ddat[,-nans]
    ##now set to binary matrix
    dfdat<-apply(ddat,1,function(x){
      sapply(unlist(x),function(y) !is.na(y) && y!='NaN')
    })
    
    return(dfdat)
    
  })
  if(length(expr.list)>1){
    comm.genes<-c()#rownames(expr.list[[1]])
    for(i in 1:length(expr.list))
      comm.genes<-union(comm.genes,rownames(expr.list[[i]]))
    
    full.dat<-do.call('cbind',lapply(expr.list,function(x){
      missing<-setdiff(comm.genes,rownames(x))
      print(length(missing))
      if(length(missing)>0)
        dat<-rbind(x[intersect(rownames(x),comm.genes),],
                   t(sapply(missing,function(y,x) rep(FALSE,ncol(x)),x)))
      else{
        dat<-x[intersect(rownames(x),comm.genes),]
      }
      colnames(dat)<-colnames(x)
      dat<-dat[comm.genes,]
      return(dat)
    }))
    
    
  }
  else{
    full.dat<-expr.list[[1]]
  }
  return(full.dat)
}

library(dplyr)

mpnst<-getDisMutationData("mpnst", "mskcc")
mpnst.df<-as.data.frame(mpnst)
mpnst.df$gene<-rownames(mpnst.df)

nf1<-filter(mpnst.df, gene == "NF1")


###no nf1 mutation data. bummer
  