library(UniprotR)
ProteinAccList<- c("P31749","Q9H161","P10275","P10415","O60885","P24385","P06850",
                   "P35222","O94907","Q86SJ6","Q9UH73","Q9UNE0","P12034",
"P21781","Q9UKV0","O43593","P05019","P01308","Q15306","P07288","Q7RTS7","O95678",
"o43790","Q9BXB1","Q8WWY8","P43657","P48449","P10636","P42898","P19838","P46531","Q13635",
"P35354","Q04206","Q01196","P04278","Q15465","Q13485","O95863", "P48436","P31213","P61812",
"P01375","Q9H3D4","Q8WVJ9","P11473","Q9GZT5","O00744")
A<- ConvertID(ProteinAccList , ID_from = "ACC+ID" , ID_to = "KEGG_ID"
          , directorypath = NULL)
B<- ConvertID(ProteinAccList , ID_from = "ACC+ID" , ID_to = "EGGNOG_ID"
              , directorypath = NULL)
C<- ConvertID(ProteinAccList , ID_from = "ACC+ID" , ID_to = "ENSEMBL_ID"
             , directorypath = NULL)
D<- ConvertID(ProteinAccList , ID_from = "ACC+ID" , ID_to = "PDB_ID"
              , directorypath = NULL)
E<- ConvertID(ProteinAccList , ID_from = "ACC+ID" , ID_to = "GENENAME"
              , directorypath = NULL)
Conversions<- cbind(A,B,C,D,E)
write.csv(Conversions, file = 'GeneIDConverted.csv')