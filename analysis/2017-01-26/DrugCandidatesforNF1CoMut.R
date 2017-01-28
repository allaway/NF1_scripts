library(synapseClient)
library(dplyr)

drugdat<-synTableQuery("SELECT * FROM syn7341038")
drugdat<-as.data.frame(drugdat@values)
drugdat.filtered<-filter(drugdat, Organism=="Homo sapiens")

topNF1comut.genie<-as.data.frame(synTableQuery("SELECT * FROM syn8073830")@values)
