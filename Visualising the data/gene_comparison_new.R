#Plotting distirbution of a gene compared to whole dataset 
library(data.table)
tpm_counts <- fread("../Data/Nonlogged_tpm_normalised_counts.csv", data.table = FALSE)
#Create mask
tpm_counts <- cbind(mask=rep(TRUE, nrow(tpm_counts)), tpm_counts)
#Find genes that have zero expression across all datasets and make mask value FALSE
RS <- rowSums(tpm_counts[,3:ncol(tpm_counts)])
tpm_counts$mask[which(RS == 0)] <- FALSE

gene_set <- tpm_counts[,2]

#Trying to get the one of all genes to be on the same scale as a single genes 
plot_dist <- function(X, genes){
  #X is a dataframe of gene counts = genes are in second column and mask in first column
  X <- X[,-1]
  #Genes = gene name/s that you want to look at the distribtuion of - maximum of 3 genes 
  #Histogram of all the data 
  all <- hist(log10(as.matrix(X[,2:ncol(X)])), breaks=50)
  #all$counts <- all$counts/ncol(X)-1 = divides by ncol 
  all$counts <- all$counts/nrow(X)
  #Find the chosen genes and log10 their counts and plot
  if(length(genes) == 3){
    #Make colours to plot
    cols <- rep(NA, length(genes)+1)
    cols[1] <- rgb(0,0,1,1/4) #Put blue in first position 
    cols[2] <- rgb(1,0,0,1/4)
    cols[3] <- rgb(1,0,1,1/4)
    cols[4] <- rgb(0,1,1,1/4)
    #Find chosen genes and plot 
    chosen <- log10(X[which(X[,1] %in% genes == TRUE),2:ncol(X)])
    chosen <- as.matrix(chosen)
    chosen_1 <- hist(chosen[1,], breaks=30)
    chosen_2 <- hist(chosen[2,], breaks=30)
    chosen_3 <- hist(chosen[3,], breaks=30)
    m <- max(max(all$counts), max(chosen_1$counts), max(chosen_2$counts), max(chosen_3$counts))
    #vals <- c(all$breaks[seq(1,length(all$breaks), 10)],8)
    plot(all, col=cols[1], main="Distribution of all genes against three genes", xlab="Log10 TPM counts", ylim=c(0,m))
    plot(chosen_1, col=cols[2], add=T)
    plot(chosen_2, col=cols[3], add=T)
    plot(chosen_3, col=cols[4], add=T)
    #axis(1, at = vals, labels=vals, las=0)
    legend("topright", c("all genes", paste(genes[1]), paste(genes[2]), paste(genes[3])), fill=c(cols[1], cols[2], cols[3], cols[4]))
    
  } else if(length(genes) == 2){
    #Make colours to plot
    cols <- rep(NA, length(genes)+1)
    cols[1] <- rgb(0,0,1,1/4)
    cols[2] <- rgb(1,0,0,1/4)
    cols[3] <- rgb(1,0,1,1/4)
    #Find chosen genes and plot 
    chosen <- log10(X[which(X[,1] %in% genes == TRUE),2:ncol(X)])
    chosen <- as.matrix(chosen)
    chosen_1 <- hist(chosen[1,], breaks=30)
    chosen_2 <- hist(chosen[2,], breaks=30)
    m <- max(max(all$counts), max(chosen_1$counts), max(chosen_2$counts))
    #vals <- c(all$breaks[seq(1,length(all$breaks), 10)],8)
    plot(all, col=cols[1], main="Distribution of all genes against two genes", xlab="Log10 TPM counts ", ylim=c(0,m))
    plot(chosen_1, col=cols[2], xlim=c(0,10), add=T)
    plot(chosen_2, col=cols[3], xlim=c(0,10), add=T)
    #axis(1, at = vals, labels=vals, las=0)
    legend("topright", c("all genes", paste(genes[1]), paste(genes[2])), fill=c(cols[1], cols[2], cols[3]))
    
  } else {
    #Make colours:
    cols <- rep(NA, length(genes)+1)
    cols[1] <- rgb(0,0,1,1/4)
    cols[2] <- rgb(1,0,0,1/4)
    #Find chosen genes and plot 
    chosen <- log10(X[which(X[,1] %in% genes == TRUE),2:ncol(X)])
    chosen <- as.matrix(chosen)
    chosen <- hist(chosen, breaks=30)
    m <- max(max(all$counts), max(chosen$counts))
    #vals <- c(all$breaks[seq(1,length(all$breaks), 10)],8)
    plot(all, col=cols[1], main="Distribution of all genes against single gene", xlab="Log10 TPM counts", ylim=c(0,m))
    plot(chosen, col=cols[2], add=T)
    #axis(1, at = vals, labels=vals, las=0)
    legend("topright", c("all genes", paste(genes)), fill=c(cols[1], cols[2]))
  }
}

#trying with 1 gene
eg <- gene_set[sample(1:length(gene_set),1)]
plot_dist(tpm_counts, eg)

#Trying with 2 genes 
eg <- gene_set[sample(1:length(gene_set),2)]
plot_dist(tpm_counts,eg)

#trying with 3 genes
eg <- gene_set[sample(1:length(gene_set),3)]
plot_dist(tpm_counts,eg)
