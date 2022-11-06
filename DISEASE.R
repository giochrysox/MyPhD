library(disgenet2r)
#To get the API key
disgenet_api_key <- get_disgenet_api_key(email ="xpysoxooy@hotmail.com",  password = "paliouri1" )
#To make the key available for all functions in the package
Sys.setenv(DISGENET_API_KEY= disgenet_api_key)

#Retrieving Gene-Disease Associations from DisGeNET
#Searching by SINGLE disease. For this example OMIM Id=109200 (AGA)
data1<- disease2gene(disease= "109200", vocabulary = "OMIM", database = "CURATED", score= c(0,1))
results1<-extract(data1)
write.csv(results1,'CAD-GENES-DISGENET.csv')

#Searching by MULTIPLE diseases. For this example: AGA and Hypotrichosis
diseasesOfInterest <- c("C4083212", "C0020678")
data2 <- disease2gene(disease = diseasesOfInterest, database = "CURATED", score =c(0,1),verbose  = TRUE )
results2<-extract(data2)
plot(data2, class = "Network",prop = 2)
write.csv(results2, 'AGA-HYporichosis-DisGeNET')

#Searching by multiple genes
myListOfGenes1 <- c("AKT1","ALX4","AR","BCL2","BRD4","CCND1","CRH","CTNNB1","DKK1","DSG4",
                   "EBF1","EDAR","FGF5","FGF7","HDAC9","HR","IGF1","INS","IRF4","KLK3",
                   "KRT74","KRT75","KRT86","LGR4","LIPH","LPAR6","LSS","MAPT","MTHFR","NFKB1",
"NOTCH1","PTCH1","PTGS2","RELA","RUNX1","SHBG","SHH","SMAD4","SNAI1","SOX9","SRD5A2",
"TGFB2","TNF","TP63","TWIST2","VDR","WNT1OA","WNT10B")                 
data3 <- gene2disease(gene= myListOfGenes1, score =c(0.3, 1),verbose  = TRUE)
results3 <- extract(data3)
write.csv(results3,'MyGenes-Diseases-DISGENET.csv')



--------------------------------------------------------------------------------
library(RCy3)
cytoscapePing()
cytoscapeVersionInfo()

# we start by using STRING App to create a network with our list of genes
installApp("stringApp")
getAppStatus('STRINGapp')
gene.list<-c("AKT1","ALX4","AR","BCL2","BRD4","CCND1","CRH","CTNNB1","DKK1","DSG4",
"EBF1","EDAR","FGF5","FGF7","HDAC9","HR","IGF1","INS","IRF4","KLK3",
"KRT74","KRT75","KRT86","LGR4","LIPH","LPAR6","LSS","MAPT","MTHFR","NFKB1",
"NOTCH1","PTCH1","PTGS2","RELA","RUNX1","SHBG","SHH","SMAD4","SNAI1","SOX9","SRD5A2",
"TGFB2","TNF","TP63","TWIST2","VDR","WNT1OA","WNT10B")
gene.str<-paste(gene.list, collapse=",")

string.cmd<-paste("string protein query cutoff=0.9999 limit=48 query",
                  gene.str, sep="=")
commandsRun(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/Mynetwork.sif',type='SIF')


#Then we load all the disease networks we want to study

string.cmd='string disease query disease="Androgenic alopecia"
cutoff=0.1 species="Homo sapiens" limit=100'
commandsGET(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/AGA.sif',type='SIF')

string.cmd='string disease query disease="Hypotrichosis"
cutoff=0.1 species="Homo sapiens" limit=100'
commandsGET(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/Hypotrichosis.sif',type='SIF')

string.cmd='string disease query disease="Atrichia with papular lesions"
cutoff=0.1 species="Homo sapiens" limit=100'
commandsGET(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/Atrichia.sif',type='SIF')

string.cmd='string disease query disease="Familial woolly hair syndrome"
cutoff=0.1 species="Homo sapiens" limit=100'
commandsGET(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/Wooly hair.sif',type='SIF')

string.cmd='string disease query disease="Prostate cancer"
cutoff=0.1 species="Homo sapiens" limit=100'
commandsGET(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/Prostate cancer.sif',type='SIF')

string.cmd='string disease query disease="Type 2 Diabetes Mellitus"
cutoff=0.1 species="Homo sapiens" limit=100'
commandsGET(string.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/T2DM.sif',type='SIF')

# From WikiPathways we can see pathways relevant to the keyword 'hair follicle'
library(rWikiPathways)
installApp("WikiPathways")
getAppStatus('STRINGapp')
keyword<-"hair follicle"
gbm.pathways<-findPathwaysByText(keyword)
# The Wikipedia id for Hair follicle development: induction is WP2804
wikipathways.cmd<-'wikipathways import-as-pathway id="WP2804"'
commandsRun(wikipathways.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/WP2804.sif',type='SIF')
# The Wikipedia id for Hair follicle development: organogenesis is WP2839
wikipathways.cmd<-'wikipathways import-as-pathway id="WP2839"'
commandsRun(wikipathways.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/WP2839',type='SIF')
# The Wikipedia id for Hair follicle development: cytodifferntiation is WP2840
wikipathways.cmd<-'wikipathways import-as-pathway id="WP2840"'
commandsRun(wikipathways.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/WP2840.sif',type='SIF')
# The Wikipedia id for EDA signaling in hair follicle development is WP3930
wikipathways.cmd<-'wikipathways import-as-pathway id="WP3930"'
commandsRun(wikipathways.cmd)
exportNetwork(filename='C:/Users/xpyso/Desktop/TS15/WP3930.sif',type='SIF')

