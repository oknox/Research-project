#Pergene ranking comparing with 2:11 genes
source("ranking_FUNCTIONS_TPM.R")
rm(RS, tpm_counts_nomask, tpm_counts) #Remove unnecessary objects

#Generate list of random numbers but using the same ones
ind <- list()
ind[[1]] <- sample(1:nrow(pergene$reclassified), 2)
for(i in 2:10){ #Use same random numbers 
  ind[[i]] <- c(ind[[i-1]], sample(1:nrow(pergene$reclassified), 1))
}

#Load in functions for multiple query of genes 
source("mulq_genes_FUNCTIONS.R") 
#How many want to flip by:
s <- seq(1, ncol(pergene$reclassified), 20)

library(foreach)
library(doParallel)
mycluster <- makeCluster(24)
registerDoParallel(mycluster)

#Looping over series of genes: 
res <- foreach(number=2:11, .packages = "data.table") %dopar% {
  pg <- pergene_rank(N = number, X = pergene$reclassified, gene = ind[[number-1]], 
                     type = "null")
  #Find the initial rankings 
  init_r <- runsum_rank(df = setDF(pg$ranked_genes), q = pg$queries, 
                        col_start = 2 + number + 1)
  #mutate and find ranking
  r <- mutate_and_rank(ranked_list = setDF(pg$ranked_genes), genes = ind[[number-1]], 
                       sequence = s, X = pergene$reclassified, r = "runsum")
  rbind(init_r, r)
}

stopCluster(mycluster)

save(res, file="pergene_2to11genes_random_dopar.Rdata")