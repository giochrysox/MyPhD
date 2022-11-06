rm(list=ls())
library(httr)
library(plyr)
library (jsonlite)
library(XML)
query1 <- "http://qsardb.org/repository/service/predictor/10967/104/models/rf?CC1=CCC(C=CCC(=CCC1)C)(C)C"
out1 <- GET(url=query1)
PredictedMeltingPoint <- content(out1)
PredictedMeltingPoint 
write.csv(PredictedMeltingPoint, "MP-alpha caryophyllene.csv")


# Clear the global environment
rm(list=ls())
#load the necessary
require(RCurl)


#First we need to import the SMILES of a compound
#(in this case we used resveratrol) that 
#we want to predict its Abrahams descriptors
smiles <- "CC1(C2CC1C(=C)C(=O)C2)C"

AbrahamL <- paste("http://qsardb.org/repository/service/predictor/10967/103/models/rf?",smiles,sep="")
AbrahamV <- paste("http://qsardb.org/repository/service/predictor/10967/102/models/rf?",smiles,sep="")
AbrahamB <- paste("http://qsardb.org/repository/service/predictor/10967/101/models/rf?",smiles,sep="")
AbrahamA <- paste("http://qsardb.org/repository/service/predictor/10967/100/models/rf?",smiles,sep="")
AbrahamS <- paste("http://qsardb.org/repository/service/predictor/10967/99/models/rf?",smiles,sep="")
AbrahamE <- paste("http://qsardb.org/repository/service/predictor/10967/98/models/rf?",smiles,sep="")
Abrahamv <- paste("http://qsardb.org/repository/service/predictor/10967/97/models/rf?",smiles,sep="")
Abrahamb <- paste("http://qsardb.org/repository/service/predictor/10967/96/models/rf?",smiles,sep="")
Abrahama <- paste("http://qsardb.org/repository/service/predictor/10967/95/models/rf?",smiles,sep="")
Abrahams <- paste("http://qsardb.org/repository/service/predictor/10967/94/models/rf?",smiles,sep="")
Abrahame <- paste("http://qsardb.org/repository/service/predictor/10967/93/models/rf?",smiles,sep="")
Abrahamc <- paste("http://qsardb.org/repository/service/predictor/10967/6/models/rf?",smiles,sep="")



AbrahamDESCL <- getURI(AbrahamL)
AbrahamDESCV <- getURI(AbrahamV)
AbrahamDESCB <- getURI(AbrahamB)
AbrahamDESCA <- getURI(AbrahamA)
AbrahamDESCS <- getURI(AbrahamS)
AbrahamDESCE <- getURI(AbrahamE)
AbrahamSolvCoeffv <- getURI(Abrahamv)
AbrahamSolvCoeffb <- getURI(Abrahamb)
AbrahamSolvCoeffa <- getURI(Abrahama)
AbrahamSolvCoeffs <- getURI(Abrahams)
AbrahamSolvCoeffe <- getURI(Abrahame)
AbrahamSolvCoeffc <- getURI(Abrahamc)
write.csv(c(smiles,AbrahamDESCL,AbrahamDESCV,AbrahamDESCB,AbrahamDESCA,AbrahamDESCS,AbrahamDESCE,AbrahamSolvCoeffv,AbrahamSolvCoeffb,AbrahamSolvCoeffa,AbrahamSolvCoeffs,AbrahamSolvCoeffe,AbrahamSolvCoeffc), "AbrahamsDescriptors-pinocarvone.csv")


