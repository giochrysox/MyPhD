#Clear the Global environment 
rm(list = ls())
#Set TS7 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS7")

#We use the plantsphytochem.csv file created from searching Pubmed with Tool Service 1
PLANTS_PHYTOCHEM_REF <- read.csv('plantsphytochem.csv')

#By setting plant name you may retrieve a list of phytochemicals contained
ACACIA_CONCINNA_PHYTOCHEM <- subset.data.frame(PLANTS_PHYTOCHEM_REF, Plant == "Acacia concinna")
write.csv(ACACIA_CONCINNA_PHYTOCHEM, file = "Acacia_Concinna_Phytochemicals_Refs.csv")

#By setting a phytochemical you can check in which plants it is contained

RESVERATROL_PLANTS <- subset.data.frame(PLANTS_PHYTOCHEM_REF, Phytochemical == "resveratrol")
write.csv(RESVERATROL_PLANTS, file = "Resveratrol_Plants.csv")

#By setting a phytochemical's PubchemCID (eg Resveratrol's CID is 445154)
#you can also check in which plants it is contained
PUBCHEMCID445154_PLANTS <- subset.data.frame(PLANTS_PHYTOCHEM_REF, PubChemCID == "445154")
write.csv(PUBCHEMCID445154_PLANTS, file = "CID445154_Plants.csv")
#--------------------------------------------------------------------------------
rm(list=ls())
library(Taxonstand)
#We use a csv file with plant names (which may not be the accepted name,
#may be a synonym, an unresolved name or contain spelling errors)
TestPlants <- read.csv('C:/users/xpyso/Desktop/test plants.csv')
# This code returns the accepted names of the plants and corrects typos etc
r1 <- TPL(TestPlants$Full.Name[1:7], corr = TRUE)
str(r1)
write.csv(r1,'CorrectedPlantNames.csv')
---------------------------------------------------------------------------------
rm(list=ls())
#Loads a csv file from Dr.Duke's database containing ethnobotanical plant uses
PLANTS_ETHNOBOTANICAL_USES <- read.csv ('ETHNOBOT.csv')
ACACIA_CONCINNA_USES <- subset.data.frame(PLANTS_ETHNOBOTANICAL_USES, TAXON == "Acacia concinna")
write.csv(ACACIA_CONCINNA_USES, file = "Acacia_Concinna_Ethnobotanical_Uses.csv")
#--------------------------------------------------------------------------------
rm(list =ls())
require("ritis")
library(RJSONIO)
#To search the TSN id by using the plants common name
TSNfromCommon <- search_common('ginger')
write.csv(TSNfromCommon,'TSNIDsfromCommonName.csv')

#To search the TSN id by using the plants scientific name
TSN <- search_scientific('Jatropha curcas', from = "all", wt = "json", raw =FALSE)
TSNID <- TSN$tsn
#To get the hierarchy by using the TSN id
hierarchy_full (tsn=TSNID, wt = "json", raw =FALSE)
A <- (hierarchy_full (tsn=TSNID, wt = "json", raw =FALSE))
B <- A$taxonname
write.csv(B, "ITISTaxonomy.csv")

#To get a list of synonyms using TSN id
synonym_names(tsn=TSNID, wt = "json", raw =FALSE)
C <- synonym_names(tsn=TSNID, wt = "json", raw =FALSE)
write.csv(C, "ITISSynonyms.csv")
#------------------------------------------------------------------------------
rm(list=ls())
library(rentrez)
library(rvest)
library(dplyr)
taxid <- entrez_search(db="taxonomy", term='Curcuma longa')$ids
url <- read_html('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=136217')
a <- url %>% html_nodes("dd a") %>% html_text()
b <- as.data.frame(a)
write.csv(b, "Curcuma longa.csv")
#------------------------------------------------------------------------------
rm(list=ls())
library(wikitaxa)
#Code below retrieves plant taxonomy from Wikidata 
A<-wt_wikispecies(name = "Acacia concinna",utf8 = FALSE)
C<-A$classification
write.csv(C,"AcaciaConcinnaWikitaxonomy.csv")
#--------------------------------------------------------------------------------
# Clear all data 
rm(list=ls())
# Import the necessary libraries
library(Reol)
library(rjson)
library(jsonlite)
# We may get Encyclopedia of Life info pages and synonyms for a specific plant. 
#In this case we used Acacia concinna
Acacia_concinna <- ("https://eol.org/api/search/Acacia+concinna.json")
A <- fromJSON(Acacia_concinna)
B <- A$results
write.csv (B, file = "Acacia_concinna_plantinfoEOL.csv")
#--------------------------------------------------------------------------------

#Clears the Global environment
rm(list = ls())
#Load library taxize
library(taxize)
#Retrieve common plant names from various databases
sci2comm(scinames='Acacia Concinna', db='eol')
sci2comm(scinames='Acacia Concinna', db='itis')
sci2comm(scinames='Acacia Concinna', db='ncbi')
sci2comm('Acacia Concinna', db='worms')
#Get all common names from 4 different databases
A <- sci2comm(scinames='Acacia Concinna', db='eol')
write.csv(A, file = "AcaciaConcinnaCommonNamesEOL.csv")
B <- sci2comm(scinames='Acacia Concinna', db='itis')
write.csv(B, file = "AcaciaConcinnaCommonNamesITIS.csv")
C <- sci2comm(scinames='Acacia Concinna', db='ncbi')
write.csv(C, file = "AcaciaConcinnaCommonNamesNCBI.csv")
D <- sci2comm('Acacia Concinna', db='worms')
write.csv(D, file = "AcaciaConcinnaCommonNamesWORMS.csv")

#-------------------------------------------------------------------------------
#This part provides info on plant geography
rm(list=ls())
library(bold)
res1 <- bold_seqspec(taxon='Acacia concinna')
geo1 <- res1$country
write.csv(geo1, 'Acacia concinna_Geo.csv', row.names = FALSE)
