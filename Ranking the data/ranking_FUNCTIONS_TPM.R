#Ranking functions but for TPM normalised data 
#Script of ranking functions
#Add mask in to counts table 
library(data.table)
tpm_counts <- fread("../../../Data/Nonlogged_tpm_normalised_counts.csv", data.table = F)
tpm_counts <- cbind(mask=rep(TRUE, nrow(tpm_counts)), tpm_counts)
#Find genes that have zero expression across all datasets and make mask value FALSE
RS <- rowSums(tpm_counts[,3:ncol(tpm_counts)])
tpm_counts$mask[which(RS == 0)] <- FALSE

#Apply mask before putting into threshold functions
applymask <- function(X){
  #X = data.frame of counts that has a mask column (col1) and gene column (col2)
  t <- X[(X[,1] == TRUE),]
  return(t)
}

tpm_counts_nomask <- tpm_counts
tpm_counts <- applymask(tpm_counts)

#Function to look for specific genes
lookupgene <- function(X, genes){
  #X = data.frame of counts, mask in col1 and genes in col2
  #genes = vector of specific genes want to look at
  if(length(genes) != 0){
    X <- X[match(genes, X[,2]),]
  } else {
    print("No genes provided")
  }
  return(X)
}

#Pergenethresholds:
setthresh <- function(X){ 
  #X = data.frame of counts with mask in col1 and gene names in col2
  #Find median of evey row = for every gene
  m <- apply(X[,3:ncol(X)], 1, function(x) median(as.numeric(x)))
  #Special case if threshold is 0 
  index <- which(m == 0)
  if(length(index) == 0){ #No threshold are 0
    t <- X[,3:ncol(X)] >= m
    t[t == TRUE] <- 1
    rownames(t) <- X[,2]
  } else if(sum(m) == 0){ #If they're all 0
    t <- X[,3:ncol(X)] > m
    t[t == TRUE] <- 1
    rownames(t) <- X[,2]
  } else { #Mix of 0s and non-zeros
    t <- X[-index,3:ncol(X)] >= m[-index]
    ex <- X[index,3:ncol(X)] > m[index]
    t <- rbind(t, ex)
    t[t == TRUE] <- 1
    #Put gene names back
    rownames(t) <- c(X[-index,2], X[index,2])
  }
  return(list("thresholds"=m, "reclassified"=t))
}

system.time(pergene <- setthresh(tpm_counts)) #33s

#Whole dataset threshold
globalthresh <- function(X){ 
  #X = data.frame of counts with mask in col1 and genes in col2
  #Take global median
  v <- as.vector(as.matrix(X[,3:ncol(X)]))
  m <- median(v)
  #Values above median value
  t <- X[,3:ncol(X)] >= m
  t[t == TRUE] <- 1
  #Put gene names back in 
  rownames(t) <- X[,2]
  return(list("threshold"=m, "binary"=t))
}

#system.time(global <- globalthresh(tpm_counts)) #10s

#Hamming function
hamming <- function(query, X) {
  #query = binarised gene want to compare similarity with 4K 
  if(length(query) != ncol(X)){
    query <- query[-c(1,2)]
    query <- as.numeric(query)
  }
  #X = data.frame of thresholded (0s and 1s) expression counts = genes as rownames
  #Find out how many positions differ 
  score <- apply(X, 1, function(j) sum(query != j))
  #Rank the hamming scores
  score <- cbind("Rank"=rank(score, ties.method = "first"), score)
  score <- as.data.table(score, keep.rownames = T)
  #Add on binary data
  X <- as.data.table(X, keep.rownames = T)
  #Set the key to gene names so can match the binary data to it
  setkey(score, rn)
  score <- score[X]
  #Sort into ascending order
  score <- score[order(Rank)]
  return(setDF(score))
}

#Hamming function that doesn't return full binary sequence attached onto scores and rankings 
ham2 <- function(query, X){
  #query = binarised gene want to get hamming distances for 
  #Find out how many positions differ 
  score <- apply(X, 1, function(j) sum(query != j))
  #Rank the hamming scores
  score <- as.data.frame(cbind("Rank"=rank(score, ties.method = "first"), score))
  score <- score[order(score$Rank),]
  return(setDT(score, keep.rownames = T))
}

#Function to apply pergene threshold to multiple dataset and binarise 
multiple_set <- function(multiple_dat, X){
  #Multiple dat = data.frame of counts with genes in first column and experiment as columns 
  #X = count matrix (mask in first column, genes in 2nd col, mask applied)
  
  #Pull out counts for those genes in multiple dataset
  X <- lookupgene(X, genes=multiple_dat[,1])
  #ind <- match(genes_selected, X[,1])
  #Get threshold for those genes 
  bin_all <- setthresh(X)
  thresh <- bin_all$thresholds
  #Remove gene names from multiple query set 
  #mul_q <- mul_q[,-1]
  #Special case if threshold is 0 
  index <- which(thresh == 0)
  if(length(index) == 0){ #If none have a threshold of 0
    t <- multiple_dat[,2:ncol(multiple_dat)] >= thresh
    t[t == TRUE] <- 1
    rownames(t) <- multiple_dat[,1]
  } else if(sum(thresh) == 0){ #If all the threshold are 0
    t <- X[,3:ncol(X)] > thresh
    t[t == TRUE] <- 1
    rownames(t) <- X[,2]
  } else {
    t <- multiple_dat[-index,2:ncol(multiple_dat)] >= thresh[-index]
    ex <- multiple_dat[index,2:ncol(multiple_dat)] > thresh[index]
    t <- rbind(t, ex)
    t[t == TRUE] <- 1
    #Put gene names back
    rownames(t) <- c(multiple_dat[-index,1], multiple_dat[index,1])
  }
  return(t)
}

#Mode function - finds most common value - will assign even numbers of 1s and 0s either based on whether the column is odd or even or randomly
Mode <- function(x, method = c("random", "column")) {
  c <- x[length(x)]; x <- x[-length(x)] #Store column number and remove 
  ux <- unique(x)
  if(length(ux) ==1){ #Only 1 value
    ux
  } else if(tabulate(match(x, ux))[1] == tabulate(match(x, ux))[2]){ #Same amount of 0s and 1s
    if(method == "random"){
      sample(c(0,1),1)
    } else if(method == "column"){
      if(c %% 2 == 0){
        1
      } else {
        0
      }
    }
  } else { #How many of each present and pick one with more
    ux[which.max(tabulate(match(x, ux)))]
  }
}

#Function that calculates the majority of a multiple query
Majority <- function(mulq, method_chosen){
  #Assign column numbers 
  mulq <- rbind(mulq, 1:ncol(mulq))
  #Run mode function
  res <- apply(mulq, 2, function(x) Mode(x, method = method_chosen))
  return(res)
}
