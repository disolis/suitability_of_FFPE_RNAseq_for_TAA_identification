
source("biomart.R")

###########################################################
## Number of detected protein-coding genes and Number of ##
## detected protein-coding genes consuming 25% of        ## 
## sequencing effort per sample                          ##
###########################################################

countData<-read.table('FC/FCcountsSTAR.txt', header=TRUE, check.names=FALSE)
rownames(countData) <- countData$GeneID
countData<-countData[, -1]

paired <- read.table("FC/PairedSamples.txt",header = T)
paired <- as.character(paired$SampleName)

ids <- c()

for (i in 1:length(paired)) {
  ids <- c(ids, which(colnames(countData)==paired[i]))
}

countData <- countData[,ids]

countData[,"GeneBiotype"] <- NA

for (i in 1:nrow(countData)){
  geneid <-  rownames(countData)[i]
  genebiotype <- BM$gene_biotype[which(BM$ensembl_gene_id==geneid)]
  countData$GeneBiotype[i] <- genebiotype
}
countData <- countData[which(countData$GeneBiotype=="protein_coding"),]
countData <- countData[,-which(colnames(countData)=="GeneBiotype")]


results <- data.frame(Sample=NA,
                      Detected_protein_coding_genes=NA,
                      Detected_protein_coding_genes25=NA)

for (i in 1:ncol(countData)) {
  
  Samplename <- as.character(colnames(countData)[i])
  
  ordered <- countData[,i]
  ordered <- sort(ordered, decreasing = T)
  
  ngenes <- length(which(ordered > 0))
  
  ngenes25 <- 0
  effort25 <- 0.25*sum(ordered) 
  j <- 1
  
  while (effort25 > 0) {
    effort25 <- effort25 - ordered[j]
    ngenes25 <- ngenes25 + 1
    j <- j+1
  }
  
  newrow <- data.frame(Sample=Samplename,
                        Detected_protein_coding_genes=ngenes,
                        Detected_protein_coding_genes25=ngenes25)
  
  results <- rbind(results, newrow)
}

results <- results[-1,]  