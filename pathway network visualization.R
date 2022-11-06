#PATHWAYS NETWORK VISUALIZATION
#Load the library
library(PANEV)
#get all pathways in KEGG
panev.pathList(string = NULL)->list
#This is the file containing genes implicated in anagen 
in.file="ALL48GENES.txt"
#This is our output
out.file="all48genes-WNT-NOTCH-HH"
#Here we define the pathways we are interested in...
#WNT signaling, Notch signaling pathway and Hedgehog signaling pathway
FL= c("path:map04310","path:map04330", "path:map04340")
#... and the organism
species='hsa'
#The levels
levels=2
panev.network(in.file = in.file, out.file = out.file, species = species, FL=FL, levels = levels)