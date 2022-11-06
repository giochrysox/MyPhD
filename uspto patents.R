#Clear the Global environment 
rm(list = ls())
#Set TS4 folder on the desktop as the working directory
setwd("C:/Users/xpyso/Desktop/TS4")
require(patentsview)
require(XML)
require(rvest)
A <- "http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&p=1&u=%2Fnetahtml%2FPTO%2Fsearch-bool.html&r=0&f=S&l=50&TERM1="
#THIS WHERE YOU PUT YOUR SEARCH TERM.IN THIS EXAMPLE ALOPECIA
TERM1 <- "Androgenetic"
B <- "+&FIELD1=TI&co1=AND&TERM2="
TERM2 <-"Alopecia"
C<-"&FIELD2=TI&d=PTXT"
sdata = c(A,TERM1, B, TERM2, C)
D <- paste(sdata[1],sdata[2],sdata[3],sdata[4],sdata[5], sep ='')
D
url<-D
download.file(url, destfile = "scrapedpage.html", quiet=TRUE)
content <- read_html("scrapedpage.html")
# Read and parse HTML file
htmlTreeParse("scrapedpage.html",useInternal = TRUE)->doc.html



# Extract all the paragraphs (HTML tag is p, starting at
# the root of the document). Unlist flattens the list to
# create a character vector.
doc.text = unlist(xpathApply(doc.html, '//td', xmlValue))

# Replace all \n by spaces
doc.text = gsub('\\n', ' ', doc.text)

# Join all the elements of the character vector into a single
# character string, separated by spaces
doc.text = paste(doc.text, collapse = ' ')

write.table(doc.text, 'RetrievePatents-AndrogeneticAlopecia.txt')


