# Clear the Global environment
rm(list = ls())
#load 'bio3d' library
library("bio3d")
# Retrieve online and read pdb file of the target ZMPSTE24
pdb <- read.pdb('2YPT')
# Retrieve the aminoacid sequence of the target
targetsequence <- pdbseq(pdb)
targetsequence
# Retrieve the binding sites
bs <- binding.site(pdb)
# Browse the results (aminoacid name and position)
bs$resnames
A <- bs$resnames
write.csv(A, 'ZMPSTE24-binding sites.csv')











