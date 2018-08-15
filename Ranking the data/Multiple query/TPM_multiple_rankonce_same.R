#Script for doing various multiple queries = 0 copies to 10 copies = j=1:11 
#Ranking once 

#Changed so mutating every row - not leaving unmutated in 

source("../../ranking_FUNCTIONS_TPM.R")

#How many bits to flip by:
N <- 20
#Seq of how many positions to flip
flipby <- seq(1, ncol(pergene$reclassified), by = N)
#Store results 
tpm_multiple_once_same <- matrix(data=NA, ncol=11, nrow=length(flipby))
#Select random row
eg <- sample(1:nrow(pergene$reclassified),1)
#Pull out that row 
q <- pergene$reclassified[eg,]
#Get hamming results 
res <- hamming(q, pergene$reclassified)
print(paste("bottom score is:", res[nrow(res),3], "for row:", eg))

for(j in 1:11){
  #Make multiple query set: 
  set <- rep(q, j)
  set <- matrix(set, nrow = j, byrow = T)
  #Initialise iterations 
  it <- 1
  #Loop over flipby: 
  for(i in flipby){
    #reset query set 
    set <- rep(q, j)
    set <- matrix(set, nrow = j, byrow = T)
    #Mutate the rows in the multiple query set 
    for(k in 1:j){
      #Pick random locations to mutate based on i
      s <- sample(1:ncol(set), i)
      set[k,s] <- !set[k,s]
    }
    #Compress into one query
    final <- Majority(set, "column")
    #Where is it ranked
    differ <- sum(final != res[1,4:ncol(res)])
    #Inequality testing: 
    x1 <- differ > res$score
    #Where is x1 greater than score 
    pos <- which(x1 == TRUE)
    #Condition if changes get converted back to original sequence:
    if(length(pos) == 0){
      tpm_multiple_once_same[it,j] <- 2
      it <- it + 1
    } else {
      tpm_multiple_once_same[it,j] <- pos[length(pos)] + 1
      it <- it + 1
    }
  }
}

save(tpm_multiple_once_same, file="TPM_multiple_rankOnce.Rdata")

load("TPM_multiple_rankOnce.Rdata")
flipby <- seq(1, 4229, by = 20)
library(RColorBrewer)
col <- brewer.pal(11, "Spectral")
plot(flipby, tpm_multiple_once_same[,1], pch=20, col=col[1], xlab="Bits flipped", ylab="Rank", main="Multiple query with copies")
for(i in 2:11){
  points(flipby, tpm_multiple_once_same[,i], pch=20, col=col[i])
}
legend("bottomright", legend=c(1:11), col=col, pch=20)
