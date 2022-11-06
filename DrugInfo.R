rm(list=ls())
#To retrieve DDI data using Pubchem CID. We first need to load ddis.csv
Drug_Drug_Interactions <- read.csv ('C:/Users/xpyso/Desktop/TS10/ddis.csv')
# 57363 is Finasterides PubChem CID
Finasteride_ddis<- subset.data.frame(Drug_Drug_Interactions, DRUG1.CID == "57363") 
write.csv(Finasteride_ddis, file = "Finasteride_DrugDrugInteractions.csv")
#--------------------------------------------------------------------------------
rm(list =ls())
#To retrieve DDI data from KEGG's database using KEGG ID
library(RCurl)
library(httr)
#KEGG ID D00418 is Minoxidil
query <- 'http://rest.kegg.jp/ddi/D00418'
GET(url=query)
A <- GET(url=query)
content(A, "raw")
bin <- content(A, "raw")
writeBin(bin, "Minoxidil-DDIs.csv")
#--------------------------------------------------------------------------------
#Clear the Global Environment
rm(list=ls())
#Read pubchemcid2adr which contains ADRs for drugs from MHRA
ADRs <- read.csv ('C:/Users/xpyso/Desktop/TS10/pubchem2adr.csv')
#Subset the dataframe to include drugs causing alopecia as an ADR
ALOPECIA_AS_ADR <- subset.data.frame(ADRs, Therapeutic.response.unexpected == "Alopecia")
#Write the .csv file of drugs causing alopecia as an ADR
write.csv(ALOPECIA_AS_ADR, 'alopeciaasADR.csv')

#--------------------------------------------------------------------------------
rm(list=ls())
# To retrieve adverse reported drug events from FDA database
patient <- read.table("DEMO19Q1.txt", sep = "$", header = T, fill = T, quote = "")
drug <- read.table("DRUG19Q1.txt", sep = "$", header = T, fill = T, quote = "")
reaction <- read.table("REAC19Q1.txt", sep = "$", header = T, fill = T, quote = "")
outcomes <- read.table("OUTC19Q1.txt", sep = "$", header = T, fill = T, quote = "")

# You can find individual drugs and their reported adverse events
# by specifying their names below. Brand names can also be used
df <- drug[(grepl("PROPECIA", drug$drugname, ignore.case = T) | # drug is likely to be entered as many different brand names, use this to capture them individually
              grepl("PROPECIA", drug$drugname, ignore.case = T)) & drug$drug_seq == 1, ] # drug seq 1 == suspect drug of many possible that patient is taking

df <- merge(df, reaction, by = "primaryid") # let's merge the drug file with reactions
df <- merge(df, outcomes, by = "primaryid") # we'll bring in outcomes, too
df2 <- as.data.frame(table(df$pt, df$outc_cod)) # count the instances of reactions and their outcomes
names(df2) <- c("reaction", "outcome", "count")
df2 <- df2[order(df2$count, decreasing = T), ]
write.csv(df2, "PROPECIA_FINASTERIDE-adrs.csv")
head(df2, 20) # top 20 reactions
#---------------------------------------------------------------------------------
rm(list=ls())
library(dplyr)
# Retrieves molecules indicated in a specific disease
indications  <- file.path('C://users/xpyso/Desktop/TS10', 'indications.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
Alopecia_Indication <- subset.data.frame(indications, meddra_name=='Alopecia')
AGA_Indiaction <- subset.data.frame(indications, meddra_name=='Androgenetic alopecia')
# Retrieves side effects caused by a specific molecule and their frequency
sideeffects <- file.path('C://users/xpyso/Desktop/TS10', 'side-effects.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
Alopecia_As_Side_Effect <- subset.data.frame(sideeffects, side_effect_name=='Alopecia')
AGA_As_Side_Effect <- subset.data.frame(sideeffects, side_effect_name=='Androgenetic alopecia')
sideterms <- file.path('C://users/xpyso/Desktop/TS10', 'side-effect-terms.tsv') %>% read.delim(na.strings='',header = TRUE,stringsAsFactors=FALSE)
sidefreq <- file.path('C://users/xpyso/Desktop/TS10', 'meddra_freq.tsv.gz') %>% read.delim(na.strings='',header = FALSE,stringsAsFactors=FALSE)
AlopeciaSE_Frequency <- subset.data.frame(sidefreq, V3=='C0002170')
AGASE_Frequency <- subset.data.frame(sidefreq, V3=='C0162311')
#---------------------------------------------------------------------------------
# Clear all data 
rm(list=ls())
#Read .csv file containing all drugs currently used in Dermatology
DERMATOLOGY_DRUGS <- read.csv ('C:/Users/xpyso/Desktop/TS10/DERMATOLOGY DRUGS.csv')
#Subset those used in alopecia
ALOPECIA_DRUGS <- subset.data.frame(DERMATOLOGY_DRUGS, Condition.Treated == "Alopecia")
write.csv(ALOPECIA_DRUGS, file = "CurrentlyUsedDrugs_Alopecia.csv")





