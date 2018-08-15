#Analysing the results of doing 10 reps of all the test lists 

#Plotting the single replicates for immunity as that is the largest set 
library(ggplot2)
library(tidyr)
library(ggthemes)
library(reshape2)

#Plotting every list as single 
files <- list.files(pattern = ".Rdata")
plot_list <- list()
it <- 1
for(each in files){
  #Load in the data
  load(each)
  #Label each replicate
  complete[[1]]$Group <- rep(1, nrow(complete[[1]]))
  f <- complete[[1]]
  for(i in 2:length(complete)){
    complete[[i]]$Group <- rep(i, nrow(complete[[i]]))
    f <- rbind(f, complete[[i]])
  }
  #Turn into factors 
  f$Genes <- as.character(f$Genes)
  f$Group <- as.factor(f$Group)
  f$Set <- as.factor(f$Set)
  #Merge into a long format 
  df <- melt(f[,2:5], id.vars=c("Group", "Set"))
  df <- unite(df, variables, c(variable, Set), remove = T)
  bp <- ggplot(data=df) + 
    geom_boxplot(aes(x=Group, y=value, fill=variables), position=position_dodge(1)) + ggtitle(paste(unlist(strsplit(each, "_"))[1], "list (n =", dim(f)[1]/length(complete), ")"))
  plot_list[[it]] <- bp + scale_fill_brewer(palette="Paired")
  it <- it+1
}

pdf("../../../figures/results_figures/mul_q/every_list_plotsingle.pdf", 
    width=12, height=5)
for(i in 1:7){
  print(plot_list[[i]])
}
dev.off()


#Plot boxplot for each list in one plot = make into one big dataframe - need column for the list it came from 
nlist <- c()
load("brain_10rep.Rdata")
complete[[1]]$Group <- rep("brain", nrow(complete[[1]]))
f <- complete[[1]]
for(i in 2:10){
  complete[[i]]$Group <- rep("brain", nrow(complete[[i]]))
  f <- rbind(f, complete[[i]])
}
nlist[1] <- dim(f)[1]/10
#Add on the other lists 
files <- list.files(pattern = ".Rdata")
files <- files[-1]
names <- substr(files, 1, nchar(files)-12)
it <- 1
for(each in files){
  load(each)
  for(i in 1:length(complete)){
    complete[[i]]$Group <- rep(names[it], nrow(complete[[i]]))
    f <- rbind(f, complete[[i]])
  }
  nlist <- c(nlist, dim(f[which(f$Group == names[it]),])[1]/10)
  it <- it + 1
}
#Turn into factors 
f$Genes <- as.character(f$Genes)
f$Group <- as.factor(f$Group)
f$Set <- as.factor(f$Set)

#Change the order of the groups from small to large
f$Group <- factor(f$Group, levels = c("brain", "hindgut", "head", "class3c",
                                      "class3a", "class3b", "imm"))

#Merge and plot 
df <- melt(f[,2:5], id.vars=c("Group", "Set"))
df <- unite(df, variables, c(variable, Set), remove = T)
pdf("../../../figures/results_figures/mul_q/All_lists.pdf", 
    width=12, height=5)
bp <- ggplot(data=df) + 
  geom_boxplot(aes(x=Group, y=value, fill=variables)) + ggtitle("All lists")
bp + scale_fill_brewer(palette="Paired") + scale_x_discrete(name="List", labels= c("brain\n(n=35)", "hindgut\n(n=39)", "head\n(n=48)", "class3c\n(n=182)", "class3a\n(n=244)", "class3b\n(n=330)", "immunity\n(n=358)")) + labs(y= "Rank")
dev.off()



#Plot scatterplot of everything in a 90% subset for rank consensus and rank pergene 
every_90_con <- f[which(f$Set == 90), "Rank_consensus"]
every_90_pg <- f[which(f$Set == 90), "Rank_pergene"]

plot(every_90_con, every_90_pg, main="Genes included in the 90% subset", xlim=c(min(min(every_90_con), min(every_90_pg)), max(max(every_90_con), max(every_90_pg))), pch=".", xlab="Rank consensus", ylab = "Rank per gene")

every_10_con <- f[which(f$Set == 10), "Rank_consensus"]
every_10_pg <- f[which(f$Set == 10), "Rank_pergene"]
plot(every_10_con, every_10_pg, main="Genes included in the 10% subset", xlim=c(min(min(every_10_con), min(every_10_pg)), max(max(every_10_con), max(every_10_pg))), pch=".", xlab="Rank consensus", ylab = "Rank per gene")

#With marginals 
library(ggExtra)

set90 <- as.data.frame(every_90_con)
colnames(set90)[1] <- "consensus_90"
set90$pergene_90 <- every_90_pg

p <- ggplot(set90, aes_string('consensus_90', 'pergene_90')) +
  geom_point(size=0.1) + theme_bw(15) + ggtitle("Genes included in 90 subset (All lists)") + 
  scale_x_continuous(limits = c(0, max(set90)))
ggExtra::ggMarginal(
  p,
  type = 'density',
  margins = 'both',
  size = 4.5,
  colour = 'black'
)

set10 <- as.data.frame(every_10_con)
colnames(set10)[1] <- "consensus_10"
set10$pergene_10 <- every_10_pg

p <- ggplot(set10, aes_string('consensus_10', 'pergene_10')) +
  geom_point(size=0.1) + theme_bw(15) + ggtitle("Genes included in 10 subset (All lists)") + 
  scale_x_continuous(limits = c(0, max(set10)))
ggExtra::ggMarginal(
  p,
  type = 'density',
  margins = 'both',
  size = 4.5,
  colour = 'black'
)


#Plotting the median, 25% quantile, 75% quantile and their standard deviations for each group 
sum.stat <- function(ranking_method, splitset){
  files <- list.files(pattern = ".Rdata")
  names <- substr(files, 1, nchar(files)-12)
  #store results 
  con <- matrix(data=NA, nrow=10,ncol=3)
  colnames(con) <- c("m", "q1", "q3")
  #Store final results 
  tab <- matrix(data=NA, nrow=7, ncol=6)
  colnames(tab) <- c("m", "m.sd", "q1", "q1.sd", "q3", "q3.sd")
  rownames(tab) <- names
  for(each in files){
    #Load in the data
    load(each)
    ind <- which(each == files)
    
    con[,1] <- unlist(lapply(complete, function(x) 
      quantile(x[which(x$Set == splitset),ranking_method], probs = .5)))
    con[,2] <- unlist(lapply(complete, function(x) 
      quantile(x[which(x$Set == splitset),ranking_method], probs = .25)))
    con[,3] <- unlist(lapply(complete, function(x) 
      quantile(x[which(x$Set == splitset),ranking_method], probs = .75)))
    
    it <- 1
    for(i in 1:3){
      tab[ind, it] <- mean(con[,i])
      tab[ind, it+1] <- sd(con[,i])
      it <- it+2
    }
  }
  return(tab)
}

con90 <- sum.stat(2, 90)
con10 <- sum.stat(2, 10)
pg90 <- sum.stat(3, 90)
pg10 <- sum.stat(3, 10)

#Reorder according to size 
ord <- c("brain", "hindgut", "head", "class3c", "class3a", "class3b", "imm")
ord <- match(ord,rownames(con90))
con90 <- con90[ord,]
con10 <- con10[ord,]
pg90 <- pg90[ord,]
pg10 <- pg10[ord,]

#Plotting q1, m and q3 in seperate plots and arranged into a group 
col <- RColorBrewer::brewer.pal(12, "Paired")
col <- col[c(1,2,3,4)]
pdf("../../../figures/results_figures/mul_q/q1_m_q3_diffmethods.pdf", height=10, width=7)
par(mfrow=c(3,1), mar=c(2,4,2,1))
for(j in c(3,1,5)){
  plot(1, type="n", xlab="", ylab="Rank", xaxt="n", xlim=c(0, 18), ylim=c(0, 12000), main=colnames(con90)[j])
  axis(1, at=c(0.75,3.25,5.75,8.25,10.75,13.25,15.75), labels = c("brain", "hindgut", "head", "class3c", "class3a", "class3b", "immunity"), tick=F)
  axis(1,at=c(0,1,1.5),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  axis(1,at=c(2.5,3,4),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  axis(1,at=c(5,6,6.5),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  axis(1,at=c(7.5,8,9),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  axis(1,at=c(10,11,11.5),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  axis(1,at=c(12.5,13,14),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  axis(1,at=c(15,16,16.5),col="black",line=0.35,tick=T,labels=rep("",3),lwd=2,lwd.ticks=0)
  
  k <- 0
  l <- 1
  for(i in 1:7){
    points(k, con10[i,j], pch=20, col= col[l])
    arrows(k, con10[i,j]-con10[i,j+1], k, 
           con10[i,j]+con10[i,j+1], length=0.05, angle=90, code=3, col=col[l])
    k <- k+0.5; l <- l+1
    points(k, con90[i,j], pch=20, col= col[l])
    arrows(k, con90[i,j]-con90[i,j+1], k, 
           con90[i,j]+con90[i,j+1], length=0.05, angle=90, code=3, col=col[l])
    k <- k+0.5; l <- l+1
    points(k, pg10[i,j], pch=20, col= col[l])
    arrows(k, pg10[i,j]-pg10[i,j+1], k, 
           pg10[i,j]+pg10[i,j+1], length=0.05, angle=90, code=3, col=col[l])
    k <- k+0.5; l <- l+1
    points(k, pg90[i,j], pch=20, col= col[l])
    arrows(k, pg90[i,j]-pg90[i,j+1], k, 
           pg90[i,j]+pg90[i,j+1], length=0.05, angle=90, code=3, col=col[l])
    k <- k+1; l <- 1
  }
  if(j == 5){
    legend("bottomright", legend = c("consensus 10", "consensus 90", "per-gene 10", "per-gene 90"), col = col, pch=20)
  } else {
    legend("topright", legend = c("consensus 10", "consensus 90", "per-gene 10", "per-gene 90"), col = col, pch=20)
  }
}
dev.off()

