#Script to make function so that when when given gene or genes it finds columns with the highest and lowest 10% and wordclouds the experiments

library(data.table)
#Read in TPM counts 
percent10 <- function(genes){
  #gene = can be a single FBgn gene name or multiple 
  tpm_counts <- fread("Nonlogged_tpm_normalised_counts.csv", data.table = F)
  tpm_counts <- cbind(mask=rep(TRUE, nrow(tpm_counts)), tpm_counts)
  #Load in metadata 
  meta <- metadata <- read.delim("20180501_metadata.tsv", stringsAsFactors = F)
  
  #Find where gene names are located 
  ind <- which(tpm_counts[,2] %in% genes == TRUE)
  
  if(length(genes) == 1){
    q <- tpm_counts[ind,-c(1,2)]
    #Find top 10% and bottom 10% 
    percent <- q[{x<-rank(q)/length(q);x<0.1 | x>=0.9}]
    #Find top and bottom 10% in the metadata 
  } else {
    q <- tpm_counts[ind,-c(1,2)]
    #Sum the columns
    qcol <- colSums(q)
    #Find top 10% and bottom 10% 
    percent <- qcol[{x<-rank(qcol)/length(qcol);x<0.1 | x>=0.9}]
  }
  #"author_source_name",
  nam <- c("author_tissue", "author_strain", "author_sex",  "author_genotype", "author_dev_stage")
  sub <- meta[which(meta[,1] %in% names(percent) == TRUE),which(colnames(meta) %in% nam == TRUE)]
  
  #Wordcloud of subset of metadata
  library(tm)
  library(SnowballC)
  library(wordcloud)
  library(RColorBrewer)
  
  v <- as.vector(as.matrix(sub))
  m <- v[v != "" & !is.na(v) & v != "Not Applicable" & v != "not applicable" & v != 'Not applicable']
  mdat <- Corpus(VectorSource(m))
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  mdat <- tm_map(mdat, toSpace, "/")
  mdat <- tm_map(mdat, toSpace, "@")
  mdat <- tm_map(mdat, toSpace, "\\|")
  mdat <- tm_map(mdat, toSpace, "%")
  mdat <- tm_map(mdat, content_transformer(tolower))
  #mdat <- tm_map(mdat, removeNumbers)
  mdat <- tm_map(mdat, removeWords, stopwords("english"))
  # Remove your own stop word specify your stopwords as a character vector
  mdat <- tm_map(mdat, removePunctuation)
  mdat <- tm_map(mdat, stripWhitespace)
  dtm <- TermDocumentMatrix(mdat)
  mat <- as.matrix(dtm)
  vec <- sort(rowSums(mat), decreasing = T)
  d <- data.frame(word = names(vec), freq = vec)
  head(d, 10)
  #generate the word cloud
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words = 200, random.order = F, rot.per = 0.35,
            colors = brewer.pal(8, "Dark2"))
}

x <- sample(1:nrow(tpm_counts), 1)
g <- tpm_counts[x,"FBgn"]

percent10(genes = g)

#Find genes that come ranked together 
#global <- globalthresh(tpm_counts)
x <- 14332
q <- pergene$reclassified[x,]
res <- hamming(q, pergene$reclassified)
g <- res[1:4,1]
percent10(genes=g)
