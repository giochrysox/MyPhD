#Clear the Global environment 
rm(list = ls())
#Set TS1 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS1")
#Load the libraries
library(RISmed)
library(rentrez)
library(XML)
#Manually insert your query 
search_topic <- 'androgenetic alopecia'
# Here, we may modify our search(e.g max 1000 articles from 1950-2021)
search_query <- EUtilsSummary(search_topic, retmax=1000, mindate=1950,maxdate=2021)
summary(search_query)
A<-search_query@querytranslation
write.csv(A, 'Query-MeSHTerms.csv')
# see the ids of our returned query
QueryId(search_query)

# get actual data from PubMed
records<- EUtilsGet(search_query)
class(records)
# store it as data frame
pubmed_data <- data.frame('PMID'=PMID(records),'Title'=ArticleTitle(records),'Abstract'=AbstractText(records))
head(pubmed_data,1)
pubmed_data$Abstract <- as.character(pubmed_data$Abstract)
pubmed_data$Abstract <- gsub(",", " ", pubmed_data$Abstract, fixed = TRUE)
# Visualise results in R
str(pubmed_data)
# write the file in csv
write.csv(pubmed_data,file='TS1-ANDROGENETIC ALOPECIA.csv')
