source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("InterMineR")

library(InterMineR)
listMines()

im <- initInterMine(mine=listMines()["FlyMine Beta"])
im
template <- getTemplates(im)
head(template)
#No template
#template[grep("gene", template$name, ignore.case=TRUE),]

#Need constraints to use package 
#From vignette: 
constraints = setConstraints(
  paths = c("Gene"),
  operators = c("LOOKUP"),
  values = list("FBgn0000150")
)

# set InterMineR-class query
queryGeneCount = setQuery(
  select = c("Gene.primaryIdentifier",
             "Gene.uberFlyRNASeqResults.count"),
  where = constraints
)
summary(queryGeneCount)
res <- runQuery(im, queryGeneCount) #Working
