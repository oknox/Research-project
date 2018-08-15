#Plotting random genes mutated with pergene method 
#Look at results
load("pergene_2to11genes_random_dopar.Rdata")
library(RColorBrewer)
cols <- brewer.pal(11,"Spectral")
s <- c(0, seq(1, 4229, 20))

pdf("../../../figures/results_figures/mul_q/pergene_ranking_random.pdf", width=10)
par(mfrow=c(2,2), mar=c(4, 4, 2.1, 5.1), xpd=TRUE)
#Plotting each gene's rank 
for(i in 1:length(res)){
  plot(s, res[[i]][,1], pch=20, xlab="Bits flipped", ylab="Rank", 
       col=cols[1], main=paste("Number of genes =", i+1), 
       ylim=c(0,16445))
  for(j in 2:ncol(res[[i]])){
    points(s, res[[i]][,j], pch=20, col=cols[j])
  }
  legend("bottomright", legend=c(1:(i+1)), pch=20, col=cols[1:(i+1)], inset = c(-0.13,0))
}
dev.off()

av <- lapply(res, function(x) rowMeans(x))
par(mar=c(5.1, 4.1, 4.1, 5.1), xpd=TRUE)
plot(s, av[[1]], pch=20, xlab="Bits flipped", ylab="Average rank", col=cols[1],
     main="Mean across all positions for differing number of genes", ylim=c(0,16445))
for(i in 2:length(av)){
  points(s, av[[i]], pch=20, col=cols[i])
}
legend("bottomright", inset=c(-0.1,0.5), legend = 2:11, col=cols[1:10], pch=20)  
