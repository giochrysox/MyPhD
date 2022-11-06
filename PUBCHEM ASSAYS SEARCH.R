# Clear the Global Environment
rm(list=ls())
#Set TS6 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS6")
#Load libraries
library(glue)
library (XML)
library(jsonlite)

#If you know the CID of a chemical and want to find assays for it.
CID<-'2244' 
#This is the PubChem CID for aspirin

url<-glue("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{CID}/aids/JSON?aids_type=active")
CID_ASSAYS<-fromJSON(url)
CID_ASSAYS<-as.data.frame(CID_ASSAYS, row.names = TRUE)
A<-CID_ASSAYS$InformationList.Information.AID[[1]]
write.csv(A, 'CID2244_Assays.TXT')

# Full records of these assays can be retrieved
AID<-'1195,1811,92967'
url2<-glue("http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{AID}/concise/JSON")
url2
ASSAYS_SUMMARY<-fromJSON(url2)
Columns <- c(ASSAYS_SUMMARY[["Table"]][["Columns"]][["Column"]])
Rows <- c(ASSAYS_SUMMARY[["Table"]][["Row"]][["Cell"]])
df <- data.frame(Rows, Columns)
write.csv(df, 'CID2244_FullAssays.csv')




