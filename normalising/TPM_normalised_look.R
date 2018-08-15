#Looking at TPM normalised counts 
setwd("Library/Mobile Documents/com~apple~CloudDocs/CAMBRIDGE/Research project/")
library(data.table)
norm_counts <- fread("Transcripts/Nonlogged_tpm_normalised_counts.csv", data.table = F)
jcounts <- norm_counts[,2:ncol(norm_counts)]
jcounts <- as.matrix(jcounts)

#Make histogram of logged 10 normalised counts 
hist(log10(jcounts), breaks = 100, xlab="log10 TPM normalised counts", main = "Distribution of log10 TPM normalised counts")


