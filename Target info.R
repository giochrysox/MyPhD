#Clear the global environment
rm(list = ls())
#Load KEGGREST library
library (KEGGREST)
# To retrieve a KEGG gene entry,name,symbol,organism,
# database links (eg PDB),aminoacid sequence, DNA sequence
# This example is for the AR gene (KEGG Id=367)
ENTRY2 <- keggGet("hsa:367")[[1]]$ENTRY 
NAME2 <- keggGet("hsa:367")[[1]]$NAME
SYMBOL2 <- keggGet("hsa:367")[[1]]$SYMBOL
ORGANISM2 <- keggGet("hsa:367")[[1]]$ORGANISM
DBLINKS2 <- keggGet("hsa:367")[[1]]$DBLINKS
AASeq <- keggGet("hsa:367")[[1]]$AASEQ
NTSeq <- keggGet("hsa:367")[[1]]$NTSEQ
df<-data.frame(ENTRY2,NAME2,SYMBOL2,ORGANISM2)
df2<-data.frame(t(DBLINKS2))
df3<-cbind(df,df2,AASeq,NTSeq)
# The resulting .csv contains all the info above
write.csv(df3,file = "hsa367GENEINFO.csv", row.names = T)

