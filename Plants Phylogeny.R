#Clear the Global environment 
rm(list = ls())
#Set TS8 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS8")
require(bold)
bold_seqspec(taxon='Jatropha curcas') 
bold_seqspec(taxon='Jatropha curcas', format='xml') 
bold_seqspec(taxon='Jatropha curcas', response=TRUE) 
res <- bold_seqspec(taxon='Jatropha curcas', sepfasta=TRUE) 
res$fasta
write.csv(res, "Jatropha Curcas-FASTA.csv")
#By repeating the process we may retrieve FASTA files for all 69 plants
#----------------------------------------------------------------------------------
rm(list=ls())
require(DECIPHER)
#We load all.plants.fas which contains the FASTA files for all 69 plants
fas <- "C:/Users/xpyso/Desktop/TS8/allplants.fas"
library(Biostrings)
dna <- readDNAStringSet(fas)
dna
dna<-RemoveGaps(dna) #removes gaps
DNA <- AlignSeqs(dna) # align the sequences directly without translation
DNA <- AlignTranslation(dna) # align the translation then reverse translate
#write the aligned sequences to a FASTA file
writeXStringSet(DNA, file="C:/users/xpyso/Desktop/TS8/msaPLANTS.fasta")

BrowseSeqs(DNA, highlight=0)
DNA_adjusted <- AdjustAlignment(DNA)
DNA_staggered <- StaggerAlignment(DNA)
d <- DistanceMatrix(DNA_staggered)
tree <- IdClusters(d, cutoff=10, method="NJ", showPlot=TRUE, type="dendrogram")