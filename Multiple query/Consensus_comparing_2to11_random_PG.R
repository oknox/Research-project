#Consensus comparing with 2:11 genes 
source("ranking_FUNCTIONS_TPM.R")
#load("random_values_for_consensus.Rdata") = seeing if pergene improves the improved ranking seen with the global threshold. - want to instead have a series of random values that are the same and ony adding on one new one each time.
#Put numbers for each number of genes into a list
#ind <- list()
#for(i in 2:11){
  #ind[[i-1]] <- x[seq(1:i)]
  #x <- x[-seq(1:i)]
#}
rm(RS, tpm_counts_nomask, tpm_counts) #Remove unnecessary objects

#Generate list of random numbers put using the same ones
load("random_numbers.Rdata")

#Load in functions for multiple query of genes 
source("mulq_genes_FUNCTIONS.R") 
#How many want to flip by:
s <- seq(1, ncol(pergene$reclassified), 20)

#Looping over series of genes: 
#Make empty list to store results 
res <- list()
for(number in 2:11){
  CR <- consensus_rank(N = number, X = pergene$reclassified, method = "random", 
                       gene = ind[[number-1]], type="null")
  r <- mutate_and_rank(ranked_list = CR$Ranked_list, genes = CR$genes, sequence = s, 
                       X = pergene$reclassified)
  r <- rbind(CR$Initial_rankings, r)
  res[[number-1]] <- r
}

save(res, file="Consensus_2to11genes_random.Rdata")

#Look at results
load("Consensus_2to11genes_random.Rdata")
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")
s <- c(0, seq(1, 4229, 20))

pdf("../../../figures/results_figures/mul_q/Consensus_ranking_random.pdf", width=10)
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

#Find which gene is getting an increase in rank 
load("../../random_values_for_consensus.Rdata")
ind <- list()
for(i in 2:11){
  ind[[i-1]] <- x[seq(1:i)]
  x <- x[-seq(1:i)]
}
#Gene improving in rank = 15437 = look at distribution
g <- rownames(pergene$reclassified)[15437]
#Load function from gene_to_whole_comparison script
plot_dist(tpm_counts, g)

#Plotting average rank at however many bits flipped
#Find rowMeans for each 
av <- lapply(res, function(x) rowMeans(x))
plot(s, av[[1]], pch=20, xlab="Bits flipped", ylab="Average rank", col=cols[1],
     type="b", main="Mean across all positions for differing number of genes", ylim=c(0,16445))
for(i in 2:length(av)){
  points(s, av[[i]], pch=20, col=cols[i], type="b")
}
legend("bottomright", legend = 2:11, col=cols[1:10], pch=20)  
