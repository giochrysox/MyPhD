require(XML)
require(plyr)

file<-'http://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot?uniprot=O00744;100'
doc <- xmlParse(file,useInternalNodes = TRUE)
xL <- xmlToList(doc)
data <- ldply(xL, data.frame)
write.csv(data,'BindingDB-WNT1OB.csv')