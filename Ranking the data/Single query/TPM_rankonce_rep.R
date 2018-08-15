#Script for doing single query with TPM data, ranking once but with replicates
source("scripts/Ranking the data/ranking_FUNCTIONS_TPM.R")

#generate sequence 
s <- seq(1,ncol(pergene$reclassified), by=20)
#create matrix to store results 
tpm_once_final_rep <- matrix(data=NA, ncol=5, nrow=length(s))

for(j in 1:5){
  #Pick random row
  ind <- sample(1:nrow(pergene$reclassified),1)
  #Pull out the row
  q <- pergene$reclassified[ind,]
  #Get the ranked list compared to the chosen row
  res <- hamming(q, pergene$reclassified)
  #Find the bottom score
  print(paste("bottom score is:", res[nrow(res),3], "for row:", ind))
  it <- 1
  #Loop over s:
  for(i in s){
    #Reset chosen row
    q <- pergene$reclassified[ind,]
    #Pick random locations to flip based on i
    flip <- sample(1:length(q), i)
    #Flip them
    q[flip] <- !q[flip]
    #How many bits different is it to the original row 
    differ <- sum(q != res[1,4:ncol(res)])
    #Where should differ be ranked:
    x1 <- differ > res$score
    #Where x1 is greater than the scores
    pos <- which(x1 == TRUE)
    #Ranking assigned
    tpm_once_final_rep[it,j] <- pos[length(pos)] + 1
    it <- it + 1
  }
}

save(tpm_once_final_rep, file="TPM_rankOnce_rep.Rdata")

load("TPM_rankOnce_rep.Rdata")
#Find which genes were used 
geneschosen <- rownames(pergene$reclassified)[c(6377,963,14338, 6687, 15958)]
s <- seq(1, 4229, by=20)
library(RColorBrewer)
col <- brewer.pal(5, "Dark2")
plot(s, tpm_once_final_rep[,1], xlab="Bits flipped", ylab="Rank", pch=20, main="Single query with replicates", col=col[1])
for(i in 2:5){
  points(s, tpm_once_final_rep[,i], col=col[i], pch=20)
}
legend("bottomright", legend =  geneschosen, pch=20, col=col)
