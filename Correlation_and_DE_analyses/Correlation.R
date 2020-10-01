##CORRELATION ANALYSIS
##INTRA-GROUP VARIABILITY
##INTER-GROUP VARIABILITY

##Include packages (install if necessary)
list.of.packages <- c("ggcorrplot", "corrplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("ggcorrplot")
library("corrplot")


##DEFINE FUNCTIONS

TPM <- function (CountsVector, LengthVector) {
  LengthVector <- LengthVector/1000 #convert gene length into kb
  PMF <- sum(CountsVector/LengthVector)/1000000 #PMF - Per Million Factor
  TPM <- CountsVector/LengthVector/PMF 
  return (TPM)
  
}

IntraGroupVar <- function(targetGroup, countDataTPM, groups) {
  ids <- which(groups==targetGroup)
  countDataTPM_oneGroup<-countDataTPM[,ids]
  corr<- cor(countDataTPM_oneGroup)
  p<-ggcorrplot(corr, hc.order = TRUE, type = "lower",lab=TRUE,
                outline.col = "white",title=paste("Intra-group variability of",targetGroup,"RNA-seq",sep=" "),
                ggtheme = ggplot2::theme_gray,
                colors = c("#6D9EC1", "white", "#E46726"))
  
  plot(p)
  # ggsave(paste("Results/","Intra_variability_of_",targetGroup,"_RNAseq.png",sep=""),plot=p,
  #        width = 5.53,
  #        height = 5.53,)
}

InterGroupVar <- function(targetGroup1,targetGroup2, countDataTPM, groups) {

  ids <- which(groups==targetGroup1)
  countDataTPM1<-countDataTPM[,ids]
  ids <- which(groups==targetGroup2)
  countDataTPM2<-countDataTPM[,ids]
  
  corr<- cor(countDataTPM2,countDataTPM1)
  
  col <- colorRampPalette(c("#77AADD", "#4477AA","#FFFFFF","#BB4444", "#EE9988"))
  # png(paste("Results/","Inter_group_variability_of_",targetGroup1,"_versus_",targetGroup2,".png",sep=""), width =600, height = 400)
  corrplot(corr, method="color",tl.srt=45,col=col(200),
           addCoef.col = "black",type="upper",insig = "blank",mar=c(0,0,1,0),title=paste("Inter-group variability of",targetGroup1,"versus",targetGroup2,sep=" "))
  # dev.off()
}

countData<-read.table('INPUT/InputCountsFC_STAR.txt', header=TRUE, check.names=FALSE)
rownames(countData)<-countData$GeneID
countData<-countData[, -1] ##remove meatadata not needed for the DE analysis

##Get sample groups (FFPE_od/FFPE_control/FF_od/FF_control)
##Each sample name has the following structure: {patient id}_[FFPE/FF]_[od/control]
##The part "[FFPE/FF]_[od/contro]" implies the sample group
groups <- c()
for (i in c(2:length(colnames(countData)))) {
  conservType<-unlist(strsplit(colnames(countData)[i],"_",fixed=T))[2]
  condition<-unlist(strsplit(colnames(countData)[i],"_",fixed=T))[3]
  group<-paste(conservType,condition,sep = "_")
  groups <- c(groups,group)
}

countDataTPM <- data.frame(countData$Length)

for (i in c(2:ncol(countData))) {
  TPMcol <- TPM(countData[,i],countData$Length)
  countDataTPM <- data.frame (countDataTPM, TPMcol)
}

colnames(countDataTPM) <- colnames(countData)
rownames(countDataTPM) <- rownames(countData)
countDataTPM <- countDataTPM[,-1]



##Intra-group variability
fourConditions <- c("FFPE_od","FFPE_control","FF_od","FF_control")
for (i in 1:4) {
  targetGroup <- fourConditions[i]
  IntraGroupVar (targetGroup, countDataTPM, groups)
}

##Inter-group variability 1: FFPE_od versus FF_od
InterGroupVar ("FFPE_od","FF_od", countDataTPM, groups) 

##Inter-group variability 2: FFPE_control versus FF_control
InterGroupVar ("FFPE_control","FF_control", countDataTPM, groups) 