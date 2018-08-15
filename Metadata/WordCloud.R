#Create word cloud 
setwd("Library/Mobile Documents/com~apple~CloudDocs/CAMBRIDGE/Research project/")
# Install
install.packages("tm")  # for text mining
install.packages("SnowballC") # for text stemming
install.packages("wordcloud") # word-cloud generator 
install.packages("RColorBrewer") # color palettes
# Load
library(tm)
library(SnowballC)
library(wordcloud)
library(RColorBrewer)

#Trying with metadata
#Load metadata
metadata <- read.delim("Data/20180501_metadata.tsv")

#Test case
test <- metadata[1:20,1:20]
test <- as.vector(as.matrix(test)) #Turn into vector 
test <- test[test != "" & !is.na(test)] #remove empty spaces and NAs 

metadata_v <- as.vector(as.matrix(metadata))
m_remove <- metadata_v[metadata_v != "" & !is.na(metadata_v) & metadata_v != "Not Applicable" & metadata_v != "not applicable" & metadata_v != 'Not applicable' & metadata_v != "not Applicable"]

mdat <- Corpus(VectorSource(m_remove))
#Inspect the content of the document
inspect(mdat)

#Text transformation- performed using tm_map() = replace special characters from the text = replace "/", "@"and "|" with space. 
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
mdat <- tm_map(mdat, toSpace, "/")
mdat <- tm_map(mdat, toSpace, "@")
mdat <- tm_map(mdat, toSpace, "\\|")

#Cleaning the text
#tm_map() function removes unnecessary white space, convert to lower case and to remove common stopwords like 'the' and 'we' 
#Can remove numbers and punctuation with removeNumbers and removePunctuation arguments
#Covert to lower case
mdat <- tm_map(mdat, content_transformer(tolower))
#Remove numbers
mdat <- tm_map(mdat, removeNumbers)
#Remove english common stopwords
mdat <- tm_map(mdat, removeWords, stopwords("english"))
# Remove your own stop word
# specify your stopwords as a character vector
#mdat <- tm_map(mdat, removeWords, c("blabla1", "blabla2"))
# Remove punctuations
mdat <- tm_map(mdat, removePunctuation)
# Eliminate extra white spaces
mdat <- tm_map(mdat, stripWhitespace)
# Text stemming
check <- tm_map(mdat, stemDocument) #Cuts off quite a few words 

#Build a term-document matrix
#Document matrix is a table containing the frequency of the words. COlumn names are words and row names are documents. 
dtm <- TermDocumentMatrix(mdat)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = T)
d <- data.frame(word = names(v), freq = v)
head(d, 15)
write.csv(head(d,15), file="../wordcloud/all_metadata_freq_table.csv")

#generate the word cloud
set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words = 200, random.order = F, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"))

#Generate word cloud with text stemming 
dtm <- TermDocumentMatrix(check)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = T)
d <- data.frame(word = names(v), freq = v)
head(d, 10)

wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words = 200, random.order = F, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"))


#Word cloud having removed drosophila, melanogaster, flies,  (where is dgrp coming from)
metadata <- read.delim("Data/20180501_metadata.tsv")
metadata_v <- as.vector(as.matrix(metadata))
m_remove <- metadata_v[metadata_v != "" & !is.na(metadata_v)]
mdat <- Corpus(VectorSource(m_remove))
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
mdat <- tm_map(mdat, toSpace, "/")
mdat <- tm_map(mdat, toSpace, "@")
mdat <- tm_map(mdat, toSpace, "\\|")
mdat <- tm_map(mdat, toSpace, "%")
mdat <- tm_map(mdat, content_transformer(tolower))
#mdat <- tm_map(mdat, removeNumbers)
mdat <- tm_map(mdat, removeWords, stopwords("english"))
# Remove your own stop word specify your stopwords as a character vector
mdat <- tm_map(mdat, removeWords, c("drosophila", "melanogaster", "flies"))
mdat <- tm_map(mdat, removePunctuation)
mdat <- tm_map(mdat, stripWhitespace)
dtm <- TermDocumentMatrix(mdat)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = T)
d <- data.frame(word = names(v), freq = v)
head(d, 10)
#generate the word cloud
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words = 200, random.order = F, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"))

findFreqTerms(dtm, lowfreq = 1000)
findAssocs(dtm, terms = "female", corlimit = 0.3)
findAssocs(dtm, terms = "w1118", corlimit = 0.3)
barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
        col ="lightblue", main ="Most frequent words",
        ylab = "Word frequencies")


#Word cloud on metadata columns that only have abstract/description/title info 
metadata <- read.delim("Data/20180501_metadata.tsv")
metadata <- metadata[,c(147,149,150)]
metadata_v <- as.vector(as.matrix(metadata))
m_remove <- metadata_v[metadata_v != "" & !is.na(metadata_v)]
mdat <- Corpus(VectorSource(m_remove))
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
mdat <- tm_map(mdat, toSpace, "/")
mdat <- tm_map(mdat, toSpace, "@")
mdat <- tm_map(mdat, toSpace, "\\|")
mdat <- tm_map(mdat, toSpace, "%")
mdat <- tm_map(mdat, content_transformer(tolower))
#mdat <- tm_map(mdat, removeNumbers)
mdat <- tm_map(mdat, removeWords, stopwords("english"))
# Remove your own stop word specify your stopwords as a character vector
mdat <- tm_map(mdat, removeWords, c("drosophila", "melanogaster", "flies"))
mdat <- tm_map(mdat, removePunctuation)
mdat <- tm_map(mdat, stripWhitespace)
dtm <- TermDocumentMatrix(mdat)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = T)
d <- data.frame(word = names(v), freq = v)
head(d, 10)
#generate the word cloud
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words = 200, random.order = F, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"))

#Word cloud on just the columns investigated not abstract/description - with and without removing numbers. Could take out column on treatment (145) and sample_extraction_method (114) and organism if high in words 
cols <- c(4,6,7,9,11,25,26,28,30,35,41,44,51,61,63,66,72,71,74,76,78,82,84,89,90,
          102,103,105,106,108,113,114,117,120,141,142,143,144,145)
descriptive_cols <- c(147,149,250)
cols_minus3 <- c(9,11,25,26,28,30,35,41,44,51,61,63,66,72,71,74,76,78,82,84,89,90,
          102,103,105,106,108,113,114,117,120,141,142,143,144,145)

metadata <- read.delim("Data/20180501_metadata.tsv")
metadata <- metadata[,cols_minus3]
metadata_v <- as.vector(as.matrix(metadata))
m_remove <- metadata_v[metadata_v != "" & !is.na(metadata_v)]
mdat <- Corpus(VectorSource(m_remove))
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
mdat <- tm_map(mdat, toSpace, "/")
mdat <- tm_map(mdat, toSpace, "@")
mdat <- tm_map(mdat, toSpace, "\\|")
mdat <- tm_map(mdat, toSpace, "%")
mdat <- tm_map(mdat, content_transformer(tolower))
#mdat <- tm_map(mdat, removeNumbers)
mdat <- tm_map(mdat, removeWords, stopwords("english"))
# Remove your own stop word specify your stopwords as a character vector
#mdat <- tm_map(mdat, removeWords, c("drosophila", "melanogaster", "flies"))
mdat <- tm_map(mdat, removePunctuation)
mdat <- tm_map(mdat, stripWhitespace)
dtm <- TermDocumentMatrix(mdat)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = T)
d <- data.frame(word = names(v), freq = v)
head(d, 10)
write.csv(head(d,15), file = "../wordcloud/Selected_columns_word_freq.csv")
#generate the word cloud
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words = 200, random.order = F, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"))

#Find where 'whole' comes up in
apply(metadata, 2, function(x) grep("whole", x)) 
apply(metadata, 2, function(x) grep("instar", x)) 

findFreqTerms(dtm, lowfreq = 300)
findAssocs(dtm, terms = "whole", corlimit = 0.3)
findAssocs(dtm, terms = "larvae", corlimit = 0.3)
barplot(d[1:10,]$freq, las = 2, names.arg = d[1:10,]$word,
        col ="lightblue", main ="Most frequent words",
        ylab = "Word frequencies")


#Word cloud on specific columns 
cols <- c(4,7,25,28,30,44,51,66,76,78,82,84,108,120,142,143,144,145)
#Make a word cloud for each of the columns - see which ones to keep 
metadata <- read.delim("Data/20180501_metadata.tsv", stringsAsFactors = F)

for(i in cols){
  meta <- metadata[,i]
  meta <- as.vector(as.matrix(meta))
  meta <- meta[meta != "" & !is.na(meta) & meta != "Not Applicable" & meta != "not applicable" & meta != 'Not applicable' & meta != "not Applicable"]
  mdat <- Corpus(VectorSource(meta))
  toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
  mdat <- tm_map(mdat, toSpace, "/")
  mdat <- tm_map(mdat, toSpace, "@")
  mdat <- tm_map(mdat, toSpace, "\\|")
  mdat <- tm_map(mdat, toSpace, "%")
  mdat <- tm_map(mdat, content_transformer(tolower))
  mdat <- tm_map(mdat, removeWords, stopwords("english"))
  mdat <- tm_map(mdat, removePunctuation)
  if(i == 145){
    mdat <- tm_map(mdat, removeNumbers)
  }
  #mdat <- tm_map(mdat, stripWhitespace)
  dtm <- TermDocumentMatrix(mdat)
  m <- as.matrix(dtm)
  v <- sort(rowSums(m), decreasing = T)
  d <- data.frame(word = names(v), freq = v)
  layout(matrix(c(1, 2), nrow=2), heights=c(0.5, 4))
  par(mar=rep(0, 4))
  plot.new()
  text(x=0.5, y=0.5, colnames(metadata)[i])
  wordcloud(words = d$word, freq = d$freq, min.freq = 1,
            max.words = 200, random.order = F, rot.per = 0.35,
            colors = brewer.pal(8, "Dark2"), main=i)
}


#Identifying how many rows are filled in the columns of wordclouds included . 
metadata <- metadata[,cols]
metadata <- metadata[,c(15,18,6)]
for(i in 1:ncol(metadata)){
  meta <- metadata[,i]
  meta <- as.vector(as.matrix(meta))
  meta <- meta[meta != "" & !is.na(meta) & meta != "Not Applicable" & meta != "not applicable" & meta != 'Not applicable' & meta != "not Applicable"]
  print(paste(length(meta), ":", colnames(metadata)[i]))
}

#Trying to merge the tissue type columns 
#Pull out the tissue columns 
tissue_cols <- c(cols[2],cols[15],cols[16],cols[17])
tissue_met <- metadata[,tissue_cols]

#Trying grepl
grepl(tissue_met[1,1], tissue_met[1,])
grepl("fat body", tissue_met[1,])

#Try adist
adist(tissue_met[17,]) #Give numbers of how different it is 

#Try agrep (levenstein distance)
#New column
tissue_met$test <- ""
for(i in 1:dim(tissue_met)[1]){
  x <- agrep(tissue_met[i,1], tissue_met[,2], ignore.case = T, value = T, max.distance = 0.05, useBytes = T)
  x <- paste0(x,"")
  tissue_met$test[i] <- x
}

#Try merge
#Remove empty rows 
store <- tissue_met == ""
rowremove <- which(rowSums(store) == 4) #Removing 1214 rows 
tissue_met <- tissue_met[-rowremove,]
#Wordcloud on author_tissue column as has the most info 
#Make column all lower case
author_tissue <- tolower(tissue_met$author_tissue)
author_tissue <- gsub("antennae", "antenna", author_tissue)
author_tissue <- gsub("antennal disc", "antenna", author_tissue)
author_tissue <- gsub("eye imaginal discs", "eye disc", author_tissue)
author_tissue <- gsub("eye-antennal disc", "eye disc", author_tissue)
author_tissue <- gsub("eye-antennal discs", "eye disc", author_tissue)
author_tissue <- gsub("eye discs", "eye disc", author_tissue)
author_tissue <- gsub("eye-antennal imaginal disc", "eye disc", author_tissue)
author_tisse <- gsub("fat body|fat bodies", "fat", author_tissue)

meta <- as.vector(as.matrix(author_tissue))
meta <- meta[meta != "" & !is.na(meta) & meta != "Not Applicable" & meta != "not applicable" & meta != 'Not applicable']
mdat <- Corpus(VectorSource(meta))
toSpace <- content_transformer(function (x , pattern ) gsub(pattern, " ", x))
mdat <- tm_map(mdat, toSpace, "/")
mdat <- tm_map(mdat, toSpace, "@")
mdat <- tm_map(mdat, toSpace, "\\|")
mdat <- tm_map(mdat, toSpace, "%")
mdat <- tm_map(mdat, content_transformer(tolower))
mdat <- tm_map(mdat, removeWords, stopwords("english"))
mdat <- tm_map(mdat, removePunctuation)
if(i == 145){
  mdat <- tm_map(mdat, removeNumbers)
}
#mdat <- tm_map(mdat, stripWhitespace)
dtm <- TermDocumentMatrix(mdat)
m <- as.matrix(dtm)
v <- sort(rowSums(m), decreasing = T)
d <- data.frame(word = names(v), freq = v)
layout(matrix(c(1, 2), nrow=2), heights=c(0.5, 4))
par(mar=rep(0, 4))
plot.new()
text(x=0.5, y=0.5, "author_tissue")
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words = 200, random.order = F, rot.per = 0.35,
          colors = brewer.pal(8, "Dark2"), main=i)



