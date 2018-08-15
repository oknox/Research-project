#Accessing the flymine lists for testing the multiple query ranking 
source("ranking_FUNCTIONS_TPM.R")
source("mulq_genes_FUNCTIONS.R")

#Loading in fly atlas brain info 
setwd("Test_fly_lists/")
hindgut_top <- read.delim("fly_atlas_hindgut.tsv", stringsAsFactors = F, header = F)
colnames(hindgut_top) <- c("Secondary_identifier", "Gene_symbol", "GeneDB_Identifier", "Organism_name")

#Make function to compare methods (consensus ranking and per gene ranking)
compare_methods <- function(tlist, X){
  #tlist = dataframe of pre-made list from FlyMine containing genes  
  #X = binarised counts 
  #Setting up the test list 
  #How many rows to get 90% 
  n = round(0.9 * nrow(tlist))
  #Which rows are going to make up the 90% 
  ind <- sample(1:nrow(tlist), n)
  #Find the selected 90% genes in the pergene dataset 
  list_names <- tlist[ind,"GeneDB_Identifier"]
  findgenes <- which(rownames(X) %in% list_names == TRUE)
  #10% left out 
  list_10 <- tlist[-ind, "GeneDB_Identifier"] 
  find10 <- which(rownames(X) %in% list_10 == TRUE)
  
  #Store the results
  foundgenes <- c(rownames(X)[findgenes], rownames(X)[find10])
  fin <- data.frame(Genes=foundgenes)
  
  ####### CONSENSUS RANKING
  #Do consensus ranking 
  r <- consensus_rank(N=NA, X=X, method="random", gene=findgenes, 
                      type = "null")
  #Rankings of the 90% inputted 
  res90 <- r$Initial_rankings
  #Find where the 10% left out are 
  y <- r$Ranked_list[which(r$Ranked_list[,"rn"] %in% list_10 == TRUE), c("rn","Rank")]
  res10 <- y$Rank; names(res10) <- y$rn ##Create named vector
  #Match the order of the gene names and insert into fin
  all <- c(res90, res10)
  fin$Rank_consensus <- all[match(names(all), fin$Genes)]
  
  ####### PERGENE RANKING
  pg <- pergene_rank(N=NA, X=X, gene = findgenes, 
                     type = "null")
  #Put 10% through runsum 
  cand <- X[find10,]
  runrank <- runsum_rank(df = setDF(pg$ranked_genes), q = cand, 
                         col_start = 2+length(findgenes)+1)
  #Put 90% through runsum 
  runrank90 <- runsum_rank(df = setDF(pg$ranked_genes), q=pg$queries,
                           col_start = 2+length(findgenes)+1)
  #Match the order of the gene names and insert into fin
  all2 <- c(runrank90, runrank)
  fin$Rank_pergene <- all2[match(names(all2), fin$Genes)]
  
  #Add to dataframe which ones are 10% and which ones are 90% 
  fin$Set <- c(rep(90, length(findgenes)), rep(10, length(find10)))
  
  return(fin)
}

#Doing 20 replicates 
complete <- list()
system.time({
  for(i in 1:10){
    complete[[i]] <- compare_methods(tlist = hindgut_top, X = pergene$reclassified)
  }
})

save(complete, file="hindgut_10rep.Rdata")