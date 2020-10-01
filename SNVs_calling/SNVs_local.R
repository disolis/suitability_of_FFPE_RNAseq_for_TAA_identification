inputDir="local/TAAs/"

## Information about the number of artefacts detected in TAAs regions for FFPE/FF pairs
## Note: in the study we focuse on C>T and G>A artefacts

files <- list.files(path = inputDir, pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


subst <- c("C>T","G>A","C>A","G>T","C>G","G>C","A>C","T>G","A>G","T>C","A>T","T>A")


listOfResults <- list()

for (k in 1:length(files)) {
  
  gene <- files[k]

  # read two times the vcf file, first for the columns names, second for the data
  tmp_vcf<-readLines(paste(inputDir,files[k],sep = ""))
  tmp_vcf_data<-read.table(paste(inputDir,files[k],sep = ""), stringsAsFactors = FALSE)

  # filter for the columns names
  tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
  vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))

  names(tmp_vcf_data)<-vcf_names

  tmp_vcf_data <- tmp_vcf_data[which(tmp_vcf_data$FILTER=="PASS"),]
  
  for (i in 1:nrow(tmp_vcf_data)) {
    for (j in 10:ncol(tmp_vcf_data)) {
      cell <- tmp_vcf_data[i,j]
      Genotype <- unlist(strsplit(cell,split = ":"))[1]
      tmp_vcf_data[i,j] <- Genotype
    }
  }

  
  n <- ncol(tmp_vcf_data)-9
  result <- matrix(0, nrow = length(subst), ncol = n)
  colnames(result) <- colnames(tmp_vcf_data)[10:ncol(tmp_vcf_data)]
  rownames(result)<-subst
    
  for (i in 1:n) {
    for (j in 1:length(subst)) {
      R<- unlist(strsplit(subst[j],split = ">"))[1]
      A<- unlist(strsplit(subst[j],split = ">"))[2]
      
      l<-length(which(tmp_vcf_data$REF==R 
                      & tmp_vcf_data$ALT==A
                      & tmp_vcf_data[,i+9] != "0/0"
                      & tmp_vcf_data[,i+9] != "./."))
      result[j,i]<-l
    }
  }
  
  gene <- strsplit(unlist(gene),split = ".",fixed = T)[[1]][1]
  listOfResults[[gene]] <- result
  # Use View() to see a particular table from the list
  # write.csv(result, paste("Local_",gene,".csv",sep = ""))

}


