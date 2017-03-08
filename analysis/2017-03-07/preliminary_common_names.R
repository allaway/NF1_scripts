library(synapseClient)
library(plyr)
library(dplyr)
library(tidyr)
synapseLogin()

syns<-read.table(synGet("syn8361358")@filePath, sep = "\t", quote = "", comment.char = "")
idx<-grep("[[:digit:]]",syns$V2)
syns2<-syns[-idx,]
syns2<-syns2[!duplicated(syns2$V1),]

syns<-filter(syns, V1 %in% syns2$V1)
syns<-bind_rows(syns, syns2)
syns<-syns[!duplicated(syns$V1),]

write.table(syns, "preliminary_commonname_map.txt", sep = "\t")
synStore(File("preliminary_commonname_map.txt", parentId = "syn7287882"))
