#Generating random genes to be used for consensus and pergene methods 

ind <- list()
ind[[1]] <- sample(1:16445, 2)
for(i in 2:10){ #Use same random numbers 
  ind[[i]] <- c(ind[[i-1]], sample(1:16445, 1))
}

save(ind, file="random_numbers.Rdata")
