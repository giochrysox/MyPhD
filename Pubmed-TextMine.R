#Clear the Global environment 
rm(list = ls())
#Set TS2 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS2")
library(pubmed.mineR)
#Here we need to download a .txt file from Pubmed containing abstracts for
#the search term
abstracts<-readabs("C://users/xpyso/desktop/TS2/ALOPECIA.txt")
abstracts@Abstract
abstrac.bodoes<-abstracts@Abstract
abstrac.bodoes
SentenceToken(abstrac.bodoes[1])
#find word
word_atomizations(abstracts)
Words<-word_atomizations(abstracts)
write.csv(Words,"ALOPECIA-WORDS.csv")
#find gene
gene_atomization(abstracts)
abs_genes = gene_atomization(abstracts)
abs_genes
write.csv(abs_genes,"ALOPECIA-GENES.csv")
abstracts@PMID
pmids<-abstracts@PMID
PMIDS<-pmids
write.csv(PMIDS,"ALOPECIA-PMIDS.csv")
pubtator_out<-pubtator_function(pmids)
pubtator_out$Genes    # shows associated genes along with their NCBI gene number
write.csv(pubtator_out$Genes ,"ALOPECIA-genes.csv")
pubtator_out$Diseases # shows associated diseases along with their MESH Id
write.csv(pubtator_out$Diseases ,"ALOPECIA-DISEASES.csv")
pubtator_out$Chemicals # shows associated chemicals along with their MESH Id
write.csv(pubtator_out$Chemicals ,"ALOPECIA-CHEMICALS.csv")
pubtator_out$Species # shows associated species along with NCBI's Taxonomy Id
write.csv(pubtator_out$Species ,"ALOPECIA-SPECIES.csv")
#--------------------------------------------------------------------------------
require(RISmed)
require(qdap)
require(tm)
require(wordcloud)
require(grDevices)
require(ggplot2)
myFunc<-function(argument){
  articles1<-data.frame('Abstract'=AbstractText(fetch), 'Year'=YearPubmed(fetch))
  abstracts1<-articles1[which(articles1$Year==argument),]
  abstracts1<-data.frame(abstracts1)
  abstractsOnly<-as.character(abstracts1$Abstract)
  abstractsOnly<-paste(abstractsOnly, sep="", collapse="")
  abstractsOnly<-as.vector(abstractsOnly)
  abstractsOnly<-strip(abstractsOnly)
  stsp<-rm_stopwords(abstractsOnly, stopwords = qdapDictionaries::Top100Words)
  ord<-as.data.frame(table(stsp))
  ord<-ord[order(ord$Freq, decreasing=TRUE),]
  head(ord,100)
}
res <- EUtilsSummary("alopecia[ti]",retmax=10000,reldate=365)
fetch <- EUtilsGet(res)
textmine2021<-myFunc(2021)
wordcloud(textmine2021$stsp, textmine2021$Freq, max.words = Inf, color = "red")

#--------------------------------------------------------------------------------
library(RISmed)

#now let's look up the term alopecia
res <- EUtilsSummary('alopecia', type='esearch', db='pubmed')
summary(res)

#what are the PubMed ids for our search?
QueryId(res)

#limit by date
res2 <- EUtilsSummary('alopecia', type='esearch', db='pubmed', mindate='1950', maxdate='2019')
summary(res2)
QueryId(res2)

tally <- array()
x <- 1
for (i in 1950:2019){
  Sys.sleep(1)
  r <- EUtilsSummary('alopecia', type='esearch', db='pubmed', mindate=i, maxdate=i)
  tally[x] <- QueryCount(r)
  x <- x + 1
}

names(tally) <- 1950:2019
max(tally)

barplot(tally, las=2, ylim=c(0,1000), main="Number of PubMed articles containing the word alopecia")