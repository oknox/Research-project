#Accessing the flymine lists for testing the multiple query ranking 
source("../../ranking_FUNCTIONS_TPM.R")
source("../../Multiple query/mulq_genes_FUNCTIONS.R")

#Loading in fly atlas brain info 
setwd("../Assessment of multiple query/Test_fly_lists/")
brain_top <- read.delim("Fly_atlas_brain_top.tsv", stringsAsFactors = F, header = F)
colnames(brain_top) <- c("Secondary_identifier", "Gene_symbol", "GeneDB_Identifier", "Organism_name")

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
  for(i in 1:20){
    complete[[i]] <- compare_methods(tlist = brain_top, X = pergene$reclassified)
  }
})
save(complete, file="brain_20rep.Rdata")


############################
load("brain_20rep.Rdata")

#Look at boxplots of single ones 
# library
library(ggplot2)

#Turning into whole dataset 
complete[[1]]$Group <- rep(1, nrow(complete[[1]]))
f <- complete[[1]]
for(i in 2:length(complete)){
  complete[[i]]$Group <- rep(i, nrow(complete[[i]]))
  f <- rbind(f, complete[[i]])
}
f$Genes <- as.character(f$Genes)
f$Group <- as.factor(f$Group)
f$Set <- as.factor(f$Set)

#Boxplots of ranking consensus results 
ggplot(f, aes(x=Group, y=Rank_consensus, fill=Set)) + 
  geom_boxplot()
#Boxplots of ranking pergene results
ggplot(f, aes(x=Group, y=Rank_pergene, fill=Set)) + 
  geom_boxplot()

library(ggplot2)
library(tidyr)
library(ggthemes)

df <- melt(f[,2:5], id.vars=c("Group", "Set"))
df <- unite(df, variables, c(variable, Set), remove = T)

pdf("../../../figures/results_figures/mul_q/brain_list_check20rep.pdf", 
    width=12, height=5)
bp <- ggplot(data=df) + 
  geom_boxplot(aes(x=Group, y=value, fill=variables), position=position_dodge(1))
bp + scale_fill_brewer(palette="Paired")
dev.off()

#Look at scatterplot of rank consensus vs rank per gene
sub90 <- f[which(f$Set == 90),]
sub10 <- f[which(f$Set == 10),]
plot(sub90$Rank_consensus, sub90$Rank_pergene, pch=20, xlab="Rank consensus", ylab="Rank pergene", main="Genes included in the 90% subset", xlim=c(min(min(sub90$Rank_consensus), min(sub90$Rank_pergene)), max(max(sub90$Rank_consensus), max(sub90$Rank_pergene))))
plot(sub10$Rank_consensus, sub10$Rank_pergene, pch=20, xlab="Rank consensus", ylab="Rank pergene", main="Genes included in the 10% subset", xlim=c(min(min(sub10$Rank_consensus), min(sub10$Rank_pergene)), max(max(sub10$Rank_consensus), max(sub10$Rank_pergene))))

#Adding marginal distributions to plot
scatterhist = function(x, y, xlab="", ylab="", main=""){
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE)
  yhist = hist(y, plot=FALSE)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  plot(x,y, pch=20, xlim=c(min(min(x), min(y)), max(max(x), max(y))))
  par(mar=c(0,3,1,1))
  barplot(xhist$counts, axes=FALSE, ylim=c(0, top), space=0)
  par(mar=c(3,0,1,1))
  barplot(yhist$counts, axes=FALSE, xlim=c(0, top), space=0, horiz=TRUE)
  par(oma=c(3,3,1,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=0.35)
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=0.35)
  #mtext(main, side=3, line=-1, outer = TRUE, adj=1,
        #at = 0.6)
  title(main, line = -1, outer = TRUE)
}
scatterhist(sub90$Rank_consensus, sub90$Rank_pergene, xlab = "Rank consensus", ylab="Rank pergene", main="Genes included in the 90% subset")
scatterhist(sub10$Rank_consensus, sub10$Rank_pergene, xlab = "Rank consensus", ylab="Rank pergene", main="Genes included in the 10% subset")

#Change to add denisty line 
install.packages("ggExtra")
library(ggExtra)

set90 <- as.data.frame(sub90$Rank_consensus)
colnames(set90)[1] <- "consensus_90"
set90$pergene_90 <- sub90$Rank_pergene

p <- ggplot(set90, aes_string('consensus_90', 'pergene_90')) +
  geom_point() + theme_bw(15) + ggtitle("Genes included in 90 subset") + 
  scale_x_continuous(limits = c(0, max(set90)))

ggExtra::ggMarginal(
  p,
  type = 'density',
  margins = 'both',
  size = 4.5,
  colour = 'black'
)

set10 <- as.data.frame(sub10$Rank_consensus)
colnames(set10)[1] <- "consensus_10"
set10$pergene_10 <- sub10$Rank_pergene

p <- ggplot(set10, aes_string('consensus_10', 'pergene_10')) +
  geom_point() + theme_bw(15) + ggtitle("Genes included in 10 subset") + 
  scale_x_continuous(limits = c(0, max(set10)))

ggExtra::ggMarginal(
  p,
  type = 'density',
  margins = 'both',
  size = 4.5,
  colour = 'black'
)

