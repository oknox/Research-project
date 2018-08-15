#Consensus comparing with 2:11 genes but selecting genes that come top with the selected gene (genes that are ranked highly together)
source("../../../ranking_FUNCTIONS_TPM.R")

#Load in functions for multiple query of genes 
source("mulq_genes_FUNCTIONS.R") 
#How many want to flip by:
s <- seq(1, ncol(pergene$reclassified), 20)

#Looping over series of genes: 
#Make empty list to store results 
res <- list()
g <- sample(1:nrow(pergene$reclassified), 1)
for(number in 2:11){
  CR <- consensus_rank(N = number, X = pergene$reclassified, method = "random", 
                       gene = g, type="close")
  r <- mutate_and_rank(ranked_list = CR$Ranked_list, genes = CR$genes, sequence = s, 
                       X = pergene$reclassified)
  r <- rbind(CR$Initial_rankings, r)
  res[[number-1]] <- r
}

save(res, file="Consensus_2to11genes_close.Rdata")

load("Consensus_2to11genes_close.Rdata")
s <- c(0, seq(1, 4229, 20))
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")

pdf("../../../figures/results_figures/mul_q/Consensus_ranking_close.pdf")
par(mfrow=c(2,2))
#Plotting each gene's rank 
for(i in 1:length(res)){
  plot(s, res[[i]][,1], pch=20, xlab="Bits flipped", ylab="Rank", 
       col=cols[1], type="b", main=paste("Number of genes =", i+1), 
       ylim=c(0,16445))
  for(j in 2:ncol(res[[i]])){
    points(s, res[[i]][,j], pch=20, type="b", col=cols[j])
  }
}
dev.off()

#Plotting average rank at however many bits flipped
#Find rowMeans for each 
av <- lapply(res, function(x) rowMeans(x))
plot(s, av[[1]], pch=20, xlab="Bits flipped", ylab="Average rank", col=cols[1],
     type="b", main="Mean across all positions for differing number of genes", ylim=c(0,16445))
for(i in 2:length(av)){
  points(s, av[[i]], pch=20, col=cols[i], type="b")
}
legend("bottomright", legend = 2:11, col=cols[1:10], pch=20)  
