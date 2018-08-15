#See how well a single gene pulls out the rest of the genes associated with it 
source("ranking_FUNCTIONS_TPM.R")

setwd("Test_fly_lists/")

files <- list.files(pattern = ".tsv")
final <- list()
for(each in files){
  #Read in file
  testlist <- read.delim(each, stringsAsFactors = F, header = F)
  colnames(testlist) <- c("Secondary_identifier", "Gene_symbol", "GeneDB_Identifier", "Organism_name")
  
  #Remove any genes not present in pergene$reclassified 
  testlist <- testlist[which(testlist$GeneDB_Identifier %in% rownames(pergene$reclassified) == TRUE),]
  
  #Pick 35 random genes 
  s <- sample(1:nrow(testlist), 35)
  
  #Loop over the randomly chosen genes and find rank of each one 
  res <- matrix(data=NA, ncol=35, nrow=dim(testlist)[1]-1)
  it <- 1
  for(i in s){
    gene <- testlist[i, "GeneDB_Identifier"]
    #Find that gene in the count table
    ind <- which(rownames(pergene$reclassified) == gene)
    #Pull it out
    q <- pergene$reclassified[ind,]
    #Rank w.r.t that gene
    dist <- ham2(q, pergene$reclassified)
    #Find where the other genes are 
    tofind <- testlist[-i, "GeneDB_Identifier"]
    res[,it] <- match(tofind, dist$rn)
    it <- it + 1
  }
  final[[which(each == files)]] <- res 
}
names(final) <- files

save(final, file="Assess_single_query_testlists.Rdata")

load("Assess_single_query_testlists.Rdata")

x <- as.vector(final[[1]])
n <- rep("class3a", length(x))
df <- data.frame(value = x, list = n)

nam <- c("class3b", "class3c", "brain", "head", "hindgut", "immunity")
for(i in 2:length(final)){
  value <- as.vector(final[[i]])
  list <- rep(nam[i-1], length(value))
  t <- cbind(value,list)
  df <- rbind(df,t)
}

df$list <- as.factor(df$list)
#Reorder the factor 
df$list <- factor(df$list, levels = c("brain", "hindgut", "head", "class3c",
                                      "class3a", "class3b", "immunity"))
df$value <- as.numeric(df$value)
pdf("../../../figures/results_figures/single_q/Boxplot_testlist.pdf", width=12, height=7)
p <- ggplot(df, aes(x = list, y = value, fill=list)) + geom_boxplot() + scale_x_discrete(name="List", labels= c("brain\n(n=35)", "hindgut\n(n=39)", "head\n(n=48)", "class3c\n(n=182)", "class3a\n(n=244)", "class3b\n(n=330)", "immunity\n(n=358)"))
p + theme(legend.position="none") + labs(y = "Rank")
dev.off()

