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