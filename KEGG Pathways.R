rm(list = ls())
library(KEGGREST)
library(png)
#Retrieve the png file of a pathway

png <- keggGet("path:hsa04310", "image")
writePNG(png, 'WNT-signaling-pathway.png')
# Mark pathway by objects and color them
url <- mark.pathway.by.objects("path:hsa04310",
                               c("hsa:1499", "hsa:595","hsa:22943","hsa:80326", 
                                 "hsa:55366","hsa:7480","hsa:4089"))
if(interactive()){
  browseURL(url)
}
