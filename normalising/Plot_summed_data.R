#Plotting sum of TPM binarised counts and sum of discretised normalised counts
setwd("../../../") #Research project working directory 

source("../../../ranking_FUNCTIONS_TPM.R")
rm(pergene, RS)
system.time(global_nomask <- globalthresh(tpm_counts_nomask))

#Sum the rows of the counts without the mask
count_nomask_sum <- rowSums(as.matrix(tpm_counts_nomask[,3:ncol(tpm_counts_nomask)]))
#Sum the rows of the binarised data 
binary_nomask_sum <- rowSums(global_nomask$binary)
plot(count_nomask_sum, binary_nomask_sum, pch=20)
plot(log10(count_nomask_sum), binary_nomask_sum, pch=20)

#Sum the rows of the counts with the mask = skewed to the left as large outlier 
count_sum <- rowSums(as.matrix(tpm_counts[,3:ncol(tpm_counts)]))
binary_sum <- rowSums(global$binary)
plot(count_sum, binary_sum, pch=20, ylab="Sum of binarised counts", xlab="Sum of counts")
plot(log10(count_sum), binary_sum, pch=20, xlab="log10 sum of counts", ylab="Sum of binarised counts")


#Test out log10 the data and colsum to see the result
#Adding in +1 as should get rid of -Inf
log10counts <- log10(as.matrix(tpm_counts[,3:ncol(tpm_counts)]) + 1) #Has -Inf 
rowsum <- rowSums(log10counts) #All -Inf 
plot(rowsum, binary_sum, pch=".", xlab="log10 TPM counts", ylab="summed binary counts")
log2counts <- log2(as.matrix(tpm_counts[,3:ncol(tpm_counts)]))
rowsum2 <- rowSums(log2counts) #Still 16,000 odd set to -Inf 

#Load in log10 data from pipeline
pipelinelog10 <- fread("../../../../Log10_tpm_normalised_counts.csv", data.table = F)
pipcolsum <- colSums(pipelinelog10[,2:ncol(pipelinelog10)])
piprowsum <- rowSums(as.matrix(pipelinelog10[,2:ncol(pipelinelog10)]))
plot(piprowsum, binary_nomask_sum, pch=".", xlab="log10 TPM counts summed", ylab="Sum of binarised counts")
