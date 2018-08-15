setwd("Library/Mobile Documents/com~apple~CloudDocs/CAMBRIDGE/Research project/Data")

counts <- read.delim("20180501_agg_gene_level_cnts.tsv", row.names = 1)
counts[1:5,1:5]

colSums(counts[,-1])

metadata <- read.delim("20180501_metadata.tsv", row.names = 1)
metadata[1:5,1:5]

#Find which experiments have metadata 
filled <- rep(NA, nrow(metadata))
names(filled) <- rownames(metadata)
for(i in 1:nrow(metadata)){
  filled[i] <- 149-length(which(metadata[i,1:149] == '' | is.na(metadata[i,1:149]) == TRUE))
}

barplot(table(filled), xlab="Number of columns containing information", ylab="Number of rows")

#For each experiment - how much metadata is present 
plot(filled, ylab="Metadata present", xlab="Experiment")
hist(filled, xlab="Metadata present", labels=TRUE, main="How much metadata does each experiment have")

#For each 'columnname/category of metadata' how much data is present across the samples 
topcategories <- rep(NA, ncol(metadata))
names(topcategories) <- colnames(metadata)[1:149]
for(j in 1:ncol(metadata)){
  topcategories[j] <- 4229-length(which(metadata[,j] == '' | is.na(metadata[,j] == TRUE)))
}

at_tick <- seq_along(topcategories)-0.5
barplot(topcategories, space=0, axes=FALSE, xaxt = 'n', horiz=TRUE)
axis(side = 2, pos=-0.2)
axis(1, at=at_tick, labels=FALSE)
axis(1, at=seq(1,149, by=5)-0.5, tick=FALSE, labels=seq(from=2,to=150,by=5), cex.lab=0.3)
axis(1, at=(1:149)-0.5, tick=FALSE, labels=names(topcategories), cex.lab=0.3, las=2)
mtext("Metadata column", side=1, line=3)
mtext("Frequency", side=2, line=1)
abline(h=500, col="red")

pdf("test_plot.pdf", height=19, width=12)
par(mar=c(1,13,1,1))
barplot(topcategories, space=0, axes=FALSE, yaxt = 'n', horiz=TRUE)
axis(side = 1, pos=-0.2)
axis(side=2, labels=names(topcategories), at=1:149-0.5, las=1)
abline(v=500, col="red")
dev.off()

#Which categories have more values in them
which(topcategories >= 500)
table(metadata[,44])


#Comparing to info from InterMineR
which(rownames(counts) == "FBgn0000150")
sum(counts[55,]) #Get 38910840

Res <- transform(res, Gene.uberFlyRNASeqResults.count = as.numeric(Gene.uberFlyRNASeqResults.count))
sum(Res[,2]) #get 108509129
