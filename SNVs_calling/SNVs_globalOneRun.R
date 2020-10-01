file="globalOneRun/od.flt.vcf"
#file="globalOneRun/control.flt.vcf"


tmp_vcf<-readLines(file)
tmp_vcf_data<-read.table(file, stringsAsFactors = FALSE)


###############################################################
## Information about the number of 12 possible substitutions ##
###############################################################

subst <- c("C>T","G>A","C>A","G>T","C>G","G>C","A>C","T>G","A>G","T>C","A>T","T>A")


# filter for the columns names
tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
  
names(tmp_vcf_data)<-vcf_names
tmp_vcf_data <- tmp_vcf_data[which(tmp_vcf_data$FILTER=="PASS"),]

#exclude 0/2, 1/2 and 2/2; there are only few 
tmp_vcf_data <- tmp_vcf_data[-which(tmp_vcf_data$ALT!="A" & 
                                    tmp_vcf_data$ALT!="T" & 
                                    tmp_vcf_data$ALT!="C" & 
                                    tmp_vcf_data$ALT!="G"),]
#Get genotype
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

result ##Result


###################################################################
## Information about the allele genotypes across all the samples ##
###################################################################

genotypes<-c("0/0","0/1","1/1","./.")

result2 <- matrix(0, nrow = length(genotypes), ncol = n)
colnames(result2) <- colnames(tmp_vcf_data)[10:ncol(tmp_vcf_data)]
rownames(result2)<-genotypes

for (i in 1:n) {
  for (j in 1:length(genotypes)) {
    l<-length(which(tmp_vcf_data[,i+9] == genotypes[j]))
    result2[j,i]<-l
  }
}

View(result2) ##Result
# write.csv(result, "globalSubstitutions.csv")
# write.csv(result2, "globalGenotypes.csv")
