# To retrieve CTD data regarding chemicals and diseases
chemicalsdiseases <- read.csv ('CTD_chemicals_diseases.csv')

alopecia_chemicals <- subset.data.frame(chemicalsdiseases, DiseaseName == "Alopecia")
write.csv(alopecia_chemicals, file = "Alopecia_Chemicals.csv")


#To retrieve CTD data regarding chemicals and genes
chemicalsgenes <- read.csv('CTD_chem_gene_ixns.csv')

AR <- subset.data.frame(chemicalsgenes, GeneSymbol == "AR")
write.csv(AR, file = "AR_Chemicals.csv")

#To retrieve CTD data regarding diseases and pathways
diseasepathways <- read.csv('CTD_diseases_pathways.csv')

alopecia <- subset.data.frame(diseasepathways, X..DiseaseName == "Alopecia")
write.csv(alopecia, file = "Alopecia_Pathways.csv")


#To retrieve CTD data regarding genes and pathways
genespathways <- read.csv('CTD_genes_pathways.csv')

AR <- subset.data.frame(genespathways, X..GeneSymbol == "AR")
write.csv(AR, file = "AR_Pathways.csv")

#To retrieve CTD data regarding AGA and GO Terms
diseaseGO <- read.csv('CTD_Disease_GO.csv')
AlopeciaGO<-subset.data.frame(diseaseGO, DiseaseName == "Alopecia")
write.csv(AlopeciaGO, 'Alopecia_GOTerms.csv')

#To retrieve CTD GO Terms related to hair follicle processes
anagen<-subset.data.frame(diseaseGO, GOID == "GO:0042640")
catagen<-subset.data.frame(diseaseGO, GOID == "GO:0042637")
haircycle<-subset.data.frame(diseaseGO, GOID == "GO:0042637")
HFdevelopment<-subset.data.frame(diseaseGO, GOID == "GO:0001942")
HFmorphogenesis<-subset.data.frame(diseaseGO, GOID == "GO:001942")
GO<-rbind(anagen,catagen,haircycle,HFdevelopment,HFmorphogenesis)
write.csv(GO, 'GO-hairprocesses.csv')
