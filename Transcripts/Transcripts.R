#Looking at gene transcripts
#Load in transcript data
tran <- read.delim("Transcripts/Transcripts.tsv", header=F, stringsAsFactors = F)
colnames(tran) <- c("DB_identifier", "Sec_identifier", "GeneSymbol", "Transcript_DB_identifier", "Transcript_length", "Organism name")
#Load in data table and the count table
library(data.table)
counts <- fread("Data/20180501_agg_gene_level_cnts.tsv", data.table = F)

#Make transcript_length numeric
tran$Transcript_length <- as.numeric(tran$Transcript_length)
#Which columns have NA's 
nacols <- function(x){
  y <- sapply(x, function(xx)any(is.na(xx)))
  names(y[y])
}  
nacols(tran) #Only Transcript_length
length(which(is.na(tran$Transcript_length) == TRUE)) #263
#Remove the rows that have NA for transcript length 
tran <- tran[-which(is.na(tran$Transcript_length) == TRUE),]

#Histogram of transcript length 
nolog <- hist(tran$Transcript_length, breaks=40, xlab="Transcript length", main="Histogram of transcript length", labels=TRUE)
plot(nolog)
#Histogram of log10 transcript length with linear values plotted - so not squished to one side 
x <- hist(log10(tran$Transcript_length), breaks=40, xlab="Transcript length", main="Histogram of transcript length")
vals <- c(x$breaks[seq(from=2,to=36,by=5)],5)
hist(log10(tran$Transcript_length), breaks=40, xlab="Transcript length", main="Histogram of transcript length", xaxt="n")
axis(1, at = vals, labels=round(10^vals), las=2)

#Calculate how many transcripts per gene and plot histogram 
genes_t <- unique(tran$DB_identifier)
library(plyr)
ntranscripts <- count(tran, "DB_identifier")
which.max(ntranscripts$freq) #What's the highest frequency 
hist(ntranscripts$freq, breaks = 40, xlab="Number of transcripts", main="Histogram of number of transcripts")
#Make histogram with log 10 number of transcripts and putting axis with linear values
y <- hist(log10(ntranscripts$freq), breaks=50, xlab="Number of transcripts", main="Histogram of number of transcripts")
vals_y <- y$breaks[c(seq(1,length(y$breaks), by=10),length(y$breaks))]
hist(log10(ntranscripts$freq), breaks=50, xlab="Number of transcripts", main="Histogram of number of transcripts", xaxt="n")
axis(1,at=vals_y, labels=round(10^vals_y))

#Collect all transcripts for a gene - calculate average transcript length and then calculate log fold change
genes_t <- unique(tran$DB_identifier)
tran$log_FC <- rep(NA, nrow(tran)) #Make column to store log2-foldchange
tran$num_assign <- rep(NA, nrow(tran)) #Make column to store number for each gene 
num <- 1
system.time({
  for(i in 1:length(genes_t)){
  #Where is that gene located 
  index <- grep(genes_t[i], tran$DB_identifier)
  tran$num_assign[index] <- num #Assign number to the genes 
  num <- num + 1
  #Mean transcript length for that gene
  m <- mean(tran[index, 5]) 
  #Calculate log2-fold change
  logFC <- log2(tran[index,5]) - log2(m)
  tran$log_FC[index] <- logFC
  }
})
plot(tran$num_assign, tran$log_FC, pch=20, xlab="Gene", ylab = "Log2FC")

#Load the means in and add to the transcript table
tran <- read.csv("Transcripts/Length_foldchange_transcripts.csv", stringsAsFactors = F)
m <- read.csv("mean_gene_lengths.csv")
tran$mean <- rep(NA, nrow(tran))
for(i in tran$num_assign){
  ind <- which(tran$num_assign == i)
  tran$mean[ind] <- m[i,2]
}
#Order according to mean length 
or <- order(tran$mean, decreasing=T)
test <- tran[or,]
#Assign number based on mean length 
test$mean_assign <- rep(NA, nrow(test))
num <- 1
genes_test <- unique(test$DB_identifier)
for(j in 1:length(genes_t)){
  index <- grep(genes_test[j], test$DB_identifier)
  test$mean_assign[index] <- num #Assign number to the genes 
  num <- num + 1
}

plot(test$mean_assign, test$log_FC, pch=20, xlab="Gene (sorted by mean transcript length in decreasing order)", ylab = "Log2FC")

#Do log10 Fold-change
#Collect all transcripts for a gene - calculate average transcript length and then calculate log fold change
genes_t <- unique(tran$DB_identifier)
tran$log10_FC <- rep(NA, nrow(tran)) #Make column to store log2-foldchange
for(i in 1:length(genes_t)){
  #Where is that gene located 
  index <- grep(genes_t[i], tran$DB_identifier)
  #tran$num_assign[index] <- num #Assign number to the genes 
  #num <- num + 1
  #Mean transcript length for that gene
  m <- mean(tran[index, 5]) 
  #Calculate log2-fold change
  logFC <- log10(tran[index,5]) - log10(m)
  tran$log10_FC[index] <- logFC
}
plot(tran$num_assign, tran$log10_FC, pch=20, xlab="Gene", ylab = "Log10FC")

m <- rep(NA, length(genes_t))
#Store mean for each of the transcripts
for(i in 1:length(genes_t)){
  index <- grep(genes_t[i], tran$DB_identifier)
  #Mean transcript length for that gene
  m[i] <- mean(tran[index, 5]) 
}

names(m) <- genes_t
write.csv(tran, file="Transcripts/Length_foldchange_transcripts.csv")
write.csv(m, file="mean_gene_lengths.csv")
