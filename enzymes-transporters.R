A<-read.csv("C:/Users/xpyso/Desktop/TS16/ENZYMES-TRANSPORTERS.csv")
Pgp <- subset.data.frame(A, METABOLISM.ENZYME.AND.TRANSPORTERS=='P-gp')
CYP3A4 <- subset.data.frame(A, METABOLISM.ENZYME.AND.TRANSPORTERS=='3A4')
PubChemCID <- subset.data.frame(A, PubChem.CID=='4201')
PubChemCIDandRef1 <- subset.data.frame(A, PubChem.CID=="4201", 
select = c('PubChem.CID','METABOLISM.ENZYME.AND.TRANSPORTERS','Relation','References'))
PubChemCIDandRef2 <- subset.data.frame(A, PubChem.CID=="57363", 
select = c(,'METABOLISM.ENZYME.AND.TRANSPORTERS','Relation','References'))
PubChemCIDandRef3 <- subset.data.frame(A, PubChem.CID=="6918296", 
select = c('PubChem.CID','METABOLISM.ENZYME.AND.TRANSPORTERS','Relation','References'))
B<-rbind(PubChemCIDandRef1,PubChemCIDandRef2,PubChemCIDandRef3)
write.csv(B, 'Drugs-Enzymes-Transporters.csv')