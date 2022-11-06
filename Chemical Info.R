#Clear the Global environment 
rm(list = ls())
#Set TS9 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS9")
library(ChemmineR)
library(ChemmineOB)
#We may retrieve SDF files from PubChem using the compoounds CID 
#e.g Curcumin=969516, Resveratrol=445154
compounds <- pubchemCidToSDF(c(969516,445154))
#SDF files can also be retrieved by inputing SMILES of a compound
sdf1 = smiles2sdf(c("COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O curcumin"))
sdf2 = smiles2sdf(c("C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O resveratrol"))
write.SDF(sdf1, file = "Curcumin.sdf")
write.SDF(sdf2, file = "Resveratrol.sdf")

#For simple chemical conversions
library(webchem)
library(jsonlite)
Name2CAS <- cir_query("Curcumin", 'cas')
CAS2SMILES <- cir_query("458-37-7", 'smiles')
Name2MW <- cir_query('Curcumin', 'mw')
# multiple inputs
comp <- c('Curcumin', 'Resveratrol')
MultipleCAS <- cir_query(comp, 'stdinchi', first = TRUE)
# to get the physicochemical properties of the compounds
URL<-'http://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/969516/property/MolecularFormula,MolecularWeight,InChIKey/JSON'
Properties<-fromJSON(URL)
Properties$PropertyTable
write.csv(Properties$PropertyTable, 'Properties.csv')
#-------------------------------------------------------------------------------
#Clear the Global environment
rm(list = ls())
#Load library rcdk
library(rcdk)
#Load the molecules we want to see if they pass Lipinski's Rule
Molecules<-load.molecules("C:/Users/xpyso/Desktop/TS9/phytochemicals.sdf")
#Generate Rule of 5 for the phytochemicals
Lipinski<-eval.desc(Molecules,
                    "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor")
#Get the molecules and number of violations
write.csv(Lipinski,"Ruleof5.csv")
#-------------------------------------------------------------------------------
#Clear Global environment
rm(list = ls())
#Load the necessary libraries
library("ChemmineR")
library("fmcsR")
#Import an sdf file
sdf <- read.SDFset("C:/users/xpyso/Desktop/TS9/drug-like phytochemicals.sdf")
#Compare a compound with a batch of compounds. In this case we compare the first
# phytochemical against all phytochemicals
sdfset <- sdf
CMP1 <- fmcsBatch(sdfset[[1]], sdfset)

write.csv(CMP1, file = "phytochemical1.csv")
#-------------------------------------------------------------------------------
#Chemical IDs Translations
library(readr)
library(dplyr)
setwd("C:/Users/xpyso/Desktop/TS9/TRANSLATION SERVICE")


#Load the translation csv's
CHEMBL<-read.csv('PUBCHEM-CHEMBL.CSV')
BINDINGDB<-read.csv('PUBCHEM-BINDINGDB.CSV')
ATC<-read.csv('PUBCHEM-ATC.CSV')
CHEBI<-read.csv('PUBCHEM-CHEBI.CSV')
DRUGBANK<-read.csv('PUBCHEM-DRUGBANK.CSV')
KEGG<-read.csv('PUBCHEM-KEGG.CSV')
#translate pubchem cid to other databases

TranslateToCHEMBL <- subset.data.frame(CHEMBL, PUBHEM.CID=='4201')
TranslateToATC <- subset.data.frame(ATC, PUBCHEM.CID=='4201')
TranslateToBINDINGDB <- subset.data.frame(BINDINGDB, PUBCHEM.CID=='4201')
TranslateToCHEBI <- subset.data.frame(CHEBI, PUBCHEM.CID=='4201')
TranslateToDRUGBANK<-subset.data.frame(DRUGBANK, PUBCHEM.CID=='4201')
TranslateToKEGG<-subset.data.frame(KEGG, PUBCHEM=='4201')

x1<- merge(TranslateToCHEMBL,TranslateToATC)  
x2<- merge(x1, TranslateToBINDINGDB)
x3<- merge(x2, TranslateToCHEBI)
x4<- merge(x3, TranslateToDRUGBANK)
x5<- merge(x4, TranslateToKEGG)
setwd("C:/Users/xpyso/Desktop/TS9")
write.csv(x5, file = "TranslationResults.csv")
#-------------------------------------------------------------------------------
#Clear the Global environment
rm(list=ls())

#Load libraries RCurl and httr
library(RCurl)
library(httr)
#Compare compound C01675 against 11,000 drugs
#and return 10 drugs with in order of maximum similarity
query <- ("http://rest.genome.jp/simcomp/C01675/drug/limit=10")
GET(url=query)
A <- GET(url=query)
content(A, "raw")
bin <- content(A, "raw")
writeBin(bin, "simcomp C01675.csv")
#-------------------------------------------------------------------------------
#Clear the Global environment
rm(list=ls())
#Load library rcdk
library(rcdk)
#Load molecules in SDF format
Molecules<-load.molecules("C:/Users/xpyso/Desktop/TS9/drug-like phytochemicals.sdf")
#Generate 2D-CDK descriptors
CDK1<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge")
CDK2<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass")
CDK3<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability")
CDK4<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor")               
CDK5<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor")
CDK6<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor")            
CDK7<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor")
CDK8<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor")
CDK9<-eval.desc(Molecules,
                "org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor")             
CDK10<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor")
CDK11<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor")
CDK12<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor")
CDK13<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor")
CDK14<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor")             
CDK15<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor")
CDK16<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor")
CDK17<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor")
CDK18<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor")
CDK19<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor")
CDK20<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor")
CDK21<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor")
CDK22<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor")                
CDK23<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor")
CDK24<-eval.desc(Molecules, 
                 "org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor")
CDK25<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor")          
CDK26<-eval.desc(Molecules,   
                 "org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor")
CDK27<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor")
CDK28<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor")
CDK29<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor")
CDK30<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor")
CDK31<-eval.desc(Molecules,     
                 "org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor")
CDK32<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor")
CDK33<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor")
CDK34<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor")
CDK35<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor")
CDK36<-eval.desc(Molecules,     
                 "org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor")        
CDK37<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor")
CDK38<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor")
CDK39<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor")
CDK40<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor")
CDK41<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor")
CDK42<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor")
CDK43<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor")
CDK44<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor")
CDK45<-eval.desc(Molecules,
                 "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor")  
CDK46<-eval.desc(Molecules, 
                 "org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor")
#Write all 2D-CDK descriptors in CSV
write.csv(c(CDK1,CDK2,CDK3,CDK4,CDK5,CDK6,CDK7,CDK8,CDK9,CDK10,CDK11,CDK12,CDK13,CDK14,CDK15,CDK16,CDK17,CDK18,CDK19,CDK20,CDK21,CDK22,CDK23,CDK24,CDK25,CDK26,CDK27,CDK28,CDK29,CDK30,CDK31,CDK32,CDK33,CDK34,CDK35,CDK36,CDK37,CDK38,CDK39,CDK40,CDK41,CDK42,CDK43,CDK44,CDK45,CDK46),"cdkall2D.csv")

#Clear the Global environment
rm(list=ls())
#Load a 3D SDF file
Molecules3D<-load.molecules("C:/Users/xpyso/Desktop/phd/methods/drug-like phytochemicals3D.sdf")
#Generate 3D-CDK descriptors
CDK3D1<-eval.desc(Molecules3D,
                  "org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor")
CDK3D2<-eval.desc(Molecules3D,
                  "org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor")
CDK3D3<-eval.desc(Molecules3D,
                  "org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor")
CDK3D4<-eval.desc(Molecules3D,
                  "org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor")
CDK3D5<-eval.desc(Molecules3D,  
                  "org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor")
CDK3D6<-eval.desc(Molecules3D,  
                  "org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor")
#Write all 3D-CDK descriptors in a CSV file
write.csv(c(CDK3D1,CDK3D2,CDK3D3,CDK3D4,CDK3D5,CDK3D6),"cdkall3D.csv")
#-------------------------------------------------------------------------------------
# To retrieve similar SMARTS patterns
library(ChemmineR)
library(ChemmineOB)
#to load an SDF file
sdfset<-read.SDFset('C:/Users/xpyso/Desktop/TS9/drug-like phytochemicals.sdf')
smartsSearchOB(sdfset[1:5],"[F,Cl,Br][CX3](=[OX1])[#1,*&!$([OH1])&!$([SH1])]",uniqueMatches=FALSE)
#or if you want to input compound's SMILES
molRefs = forEachMol("SMILES","C1CCCCC1\ttest-compound-name",identity)
smartsSearch_OB(molRefs,"[F,Cl,Br][CX3](=[OX1])[#1,*&!$([OH1])&!$([SH1])]")

