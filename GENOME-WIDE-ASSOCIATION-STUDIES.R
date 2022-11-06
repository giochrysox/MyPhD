#Clear the Global environment 
rm(list = ls())
#Set TS3 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS3")

#File gwas.csv is downloaded from EBI's GWAS Catalogue (https://www.ebi.ac.uk/gwas/docs/file-downloads)

A<-read.csv("C://users/xpyso/Desktop/TS3/gwas.csv")
#and then subset to fit androgenetic alopecia
B<-subset.data.frame(A,DISEASE.TRAIT=='Male-pattern baldness')
#We can also subset according to gene
D<-subset.data.frame(A,MAPPED_GENE=='SRD5A2')
#and store the two files in .csv format
write.csv(B,"alopecia-gwas.csv")
write.csv(D, "SRD5A2-gwas.csv")