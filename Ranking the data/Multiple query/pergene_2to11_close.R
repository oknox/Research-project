#Pergene ranking comparing with 2:11 genes
source("ranking_FUNCTIONS_TPM.R")
rm(RS, tpm_counts_nomask, tpm_counts) #Remove unnecessary objects

#Load in functions for multiple query of genes 
source("mulq_genes_FUNCTIONS.R") 
#How many want to flip by:
s <- seq(1, ncol(pergene$reclassified), 20)

library(foreach)
library(doParallel)
mycluster <- makeCluster(24)
registerDoParallel(mycluster)

g <- sample(1:nrow(pergene$reclassified), 1)
#Looping over series of genes: 
res <- foreach(number=2:11, .packages = "data.table") %dopar% {
  pg <- pergene_rank(N = number, X = pergene$reclassified, gene = g, 
                     type = "close")
  #Find the initial rankings 
  init_r <- runsum_rank(df = setDF(pg$ranked_genes), q = pg$queries, 
                        col_start = 2 + number + 1)
  #mutate and find ranking
  r <- mutate_and_rank(ranked_list = setDF(pg$ranked_genes), 
                       genes = which(rownames(pergene$reclassified) %in% rownames(pg$queries) == TRUE), sequence = s, X = pergene$reclassified, r = "runsum")
  rbind(init_r, r)
}

stopCluster(mycluster)

save(res, file="pergene_2to11genes_close_dopar.Rdata")

#Look at results
load("pergene_2to11genes_close_dopar.Rdata")
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")
s <- c(0, seq(1, 4229, 20))

pdf("../../../figures/results_figures/mul_q/Pergene_ranking_close.pdf", width=10)
par(mfrow=c(2,2), mar=c(4, 4, 2.1, 5.1), xpd=TRUE)
#Plotting each gene's rank 
for(i in 1:length(res)){
  plot(s, res[[i]][,1], pch=20, xlab="Bits flipped", ylab="Rank", 
       col=cols[1], type="b", main=paste("Number of genes =", i+1), 
       ylim=c(0,16445))
  for(j in 2:ncol(res[[i]])){
    points(s, res[[i]][,j], pch=20, col=cols[j])
  }
  legend("bottomright", legend=c(1:(i+1)), pch=20, col=cols[1:(i+1)], inset = c(-0.13,0))
}
dev.off()


#Plotting average rank at however many bits flipped
#Find rowMeans for each 
av <- lapply(res, function(x) rowMeans(x))
par(mar=c(4, 4, 2.1, 5.1), xpd=TRUE)
plot(s, av[[1]], pch=20, xlab="Bits flipped", ylab="Average rank", col=cols[1],
     type="b", main="Mean across all positions for differing number of genes", ylim=c(0,16445))
for(i in 2:length(av)){
  points(s, av[[i]], pch=20, col=cols[i], type="b")
}
legend("bottomright", legend = 2:11, col=cols[1:10], pch=20,inset = c(-0.1,0))  
