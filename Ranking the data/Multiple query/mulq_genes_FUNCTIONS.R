#Script for mulq functions 

#Function which will pick N random genes - form a consensus then rank according to that and then will return the ranking
consensus_rank <- function(N, X, method, gene, type){
  #N = how many genes want to pull out
  #X = pergene$reclassified 
  #method = whether to assign 0 or 1 by column or randomly 
  #gene = gene indices
  #type = if specify 'close' will pull out genes that get ranked to the top with the chosen gene
  #Pick genes out at random 
  if(type == "close"){
    eg  <- gene
    q <- X[eg,]
    hdist <- hamming(q, X)
    #Pick how many genes you need
    name <- hdist[1:N,1]
    eg <- which(rownames(X) %in% name == TRUE)
  } else {
    eg <- gene
  }
  q <- X[eg,]
  #Compress into majority consensus
  c <- Majority(q, method)
  #Merge consensus with ranked list
  new <- rbind(X, c)
  #Get ranked list of everything relative to the consensus
  dist <- hamming(c, new)
  #Find bottom score
  #print(paste("bottom score is:", dist[nrow(dist), 3], "for rows:", eg))
  #Find where chosen genes are ranked: 
  vec <- rep(NA, length(eg))
  names(vec) <- rownames(X)[eg]
  for(i in 1:length(eg)){
    vec[i] <- dist[which(dist$rn == rownames(X)[eg[i]]), 2]
  }
  return(list("Initial_rankings"=vec, "Ranked_list" = dist, "genes" = eg))
}

#Function will pick N random genes or given a set of positions of genes - rank pergene$reclassified by each of the genes and then combine them into one by averaging over all the gene positions 
pergene_rank <- function(N, X, gene, type){
  #N = how many genes want to pull out
  #X = pergene$reclassified 
  #gene = gene index/indices
  #Type = if want genes close by or just by indices
  #Pick out gene 
  if(type == "close"){
    eg  <- gene
    q <- X[eg,]
    hdist <- hamming(q, X)
    #Pick how many genes you need
    name <- hdist[1:N,1]
    eg <- which(rownames(X) %in% name == TRUE)
  } else {
    eg <- gene
  }
  q <- X[eg,]
  #Find ranking for each gene
  eachgene <- list()
  for(i in 1:length(eg)){
    eachgene[[i]] <- ham2(q[i,], X)
  }
  
  #Find the rank of the genes in the other list 
  f <- eachgene[[1]][,-3]
  #f <- as.data.table(f, keep.rownames=T)
  for(l in 2:length(eachgene)){
    f[,paste("list",l) := match(f$rn, eachgene[[l]]$rn)]
  }
  
  #Get average rank for the ranked lists and add to list
  f[,mean := rowMeans(f[,2:(length(eachgene)+1)])]
  #Attach on binary sequence
  X <- as.data.table(X, keep.rownames = T)
  newres <- merge(x = f, y = X, by = "rn")
  #Sort by average rank
  newres <- newres[order(newres$mean),]
  
  return(list("ranked_genes"=newres, "queries"=q))
}

#Using algorithm to find best place for 
runsum_rank <- function(df, q, col_start){
  #df = sorted list of ranked genes based on mean rank with binary sequences attached
  #q = matrix of genes and their binary counts 
  rs <- rep(NA, nrow(q))
  names(rs) <- rownames(q)
  for(i in 1:nrow(q)){
    #Calculate hamming distance from gene to ranked list 
    score <- apply(df[,col_start:ncol(df)], 1, function(j) sum(q[i,] != j))
    #Do runsum - make into separate function 
    #Initialise 
    max <- 0; #Maximum difference
    i_max <- 0; #Position of maximum difference 
    runsum <- 0 #Difference at current position
    store_runsum <- rep(NA, length(2:length(score)))
    #Insert a second first row for comparison 
    #score <- c(score[1], score)
    if(score[1] == 0){
      i_max <- 1
      runsum <- 1
    }
    for(j in 2:length(score)){
      if(score[j-1] - score[j] < 1){ 
        runsum <- runsum - 1
      } else if(score[j-1] - score[j] > 1){
        runsum <- runsum +  1
        if(runsum > max){
          max <- runsum; 
          i_max <- j
        }
      }
      store_runsum[j-1] <- runsum
    }
    rs[i] <- i_max
  }
  return(rs)
}


#Function for mutating the genes and ranking based on a sorted list
mutate_and_rank <- function(ranked_list, genes, sequence, X, r){
  #ranked_list = list that are using to rank 
  #genes = index of the genes in the counts table
  #s = sequence of values to flip by 
  #X = pergene$reclassified
  #Store results 
  ranking <- matrix(data=NA, ncol=length(genes), nrow=length(sequence))
  colnames(ranking) <- rownames(X)[genes]
  #Initialise iterations 
  it <- 1
  for(i in sequence){
    #Pull out genes
    q <- X[genes,]
    
    #Runsum method
    if(r == "runsum"){
      for(j in 1:length(genes)){
        choose <- sample(1:ncol(q), i)
        q[j,choose] <- !q[j,choose]
      }
      x <- runsum_rank(df = ranked_list, q = q, col_start = 2+length(genes)+1)
      ranking[it,] <- x
    } else {
      for(j in 1:length(genes)){
        choose <- sample(1:ncol(q), i)
        q[j,choose] <- !q[j,choose]
        #Where is it ranked
        differ <- sum(q[j,] != ranked_list[1,4:ncol(ranked_list)])
        #Inequality testing: 
        x1 <- differ > ranked_list$score
        #Where is x1 greater than score 
        pos <- which(x1 == TRUE)
        ranking[it,j] <- pos[length(pos)] + 1
      }
    }
    it <- it+1
  }
  return(ranking)
}
