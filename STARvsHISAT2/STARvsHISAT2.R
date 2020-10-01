#################################################################
## Comparative analysis of two alignment tools Hisat2 and STAR ##
#################################################################


##Include packages (install if necessary)

list.of.packages <- c("plyr","RColorBrewer","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(plyr)
library(ggplot2)
library(RColorBrewer)

source("biomart.R") ##Load biomart db
biotypes <- levels(factor(BM$gene_biotype))

########################################################
## Comparative analysis of alignment: Hisat2 vs. STAR ##
########################################################

Alignment <- read.csv("INPUT/Alignment_Hisat2_vs_STAR.csv", header=T, check.names = F)


datHISAT2 <- data.frame(Sample=Alignment$`Sample name`, 
                        Alignment=Alignment$`% Proper Pairs by Hisat2`,
                        Tool=rep("Hisat2", nrow(Alignment)))
datHISAT2$Alignment <- round(datHISAT2$Alignment*100)

datSTAR <- data.frame(Sample=Alignment$`Sample name`, 
                      Alignment=Alignment$`% Proper Pairs by STAR`,
                      Tool=rep("STAR", nrow(Alignment)))
datSTAR$Alignment <- round(datSTAR$Alignment*100)

dat <- rbind(datHISAT2,datSTAR)

#Note: uncomment ggsave() if you want to save the plot
ggplot(dat, aes(x=Sample, y=Alignment, group=Tool, color=Tool)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Hisat2 vs. STAR. % of reads mapped in proper pair") +
  geom_point(aes(shape=Tool), size=2) +
  geom_line(aes(linetype=Tool), size=0.7) 
# ggsave("Alignment_Hisat2_vs_STAR.png", width = 10, height = 5)


#################################################################
## Read input counts from two alignment tools: Hisat2 and STAR ##
#################################################################

countDataSTAR<-read.table("INPUT/FCcountsSTAR.txt", header=TRUE, check.names = F)
countDataHISAT2<-read.table("INPUT/FCcountsHISAT2.txt", header=TRUE, check.names = F)

rownames(countDataSTAR) <- countDataSTAR$GeneID
countDataSTAR <- countDataSTAR[,-1]

rownames(countDataHISAT2) <- countDataHISAT2$GeneID
countDataHISAT2 <- countDataHISAT2[,-1]

rownames(BM)<-BM$ensembl_gene_id
BM <- BM[,-c(2,4)]

countDataHISAT2 <- cbind(BM,countDataHISAT2)
countDataSTAR <- cbind(BM,countDataSTAR)

countDataHISAT2 <- countDataHISAT2[,-1]
countDataSTAR <- countDataSTAR[,-1]


##################################################
## Build distribution of counts by gene biotype ##
##################################################


dfSTAR<-data.frame("Sample"=NA,
                   "geneBiotype"=NA,
                   "Percent_of_Counts"=NA)

for (i in 2:ncol(countDataSTAR)) {
  n <- sum(countDataSTAR[,i])
  for (j in 1:length(biotypes)) {
    n1<-sum(countDataSTAR[which(countDataSTAR$gene_biotype==biotypes[j]),i])
    p<-(n1/n)*100 #percentage of counts of biotypes[j]
    newrow <- data.frame(colnames(countDataSTAR[i]),biotypes[j],p)
    names(newrow) <- c("Sample","geneBiotype","Percent_of_Counts")
    dfSTAR<-rbind(dfSTAR,newrow)
  }
}
dfSTAR<-dfSTAR[-1,]

##Split data frame, apply function, and return results in a data frame. Add position in plot label_ypos
dfSTAR <- ddply(dfSTAR, "Sample",
                   transform, label_ypos=cumsum(Percent_of_Counts))



dfHISAT2 <- data.frame("Sample"=NA,
                   "geneBiotype"=NA,
                   "Percent_of_Counts"=NA)

for (i in 2:ncol(countDataHISAT2)) {
  n <- sum(countDataHISAT2[,i])
  for (j in 1:length(biotypes)) {
    n1<-sum(countDataHISAT2[which(countDataHISAT2$gene_biotype==biotypes[j]),i])
    p<-(n1/n)*100 #percentage of counts of biotypes[j]
    newrow <- data.frame(colnames(countDataHISAT2[i]),biotypes[j],p)
    names(newrow) <- c("Sample","geneBiotype","Percent_of_Counts")
    dfHISAT2<-rbind(dfHISAT2,newrow)
  }
}
dfHISAT2<-dfHISAT2[-1,]
dfHISAT2 <- ddply(dfHISAT2, "Sample",
                  transform, label_ypos=cumsum(Percent_of_Counts))



# Define the number of colors you want
colourCount = 20
getPalette = colorRampPalette(brewer.pal(9, "Set1"))


##STAR
##Note: uncomment ggsave() if you want to save plot

ggplot(data=dfSTAR, aes(x=Sample, y=Percent_of_Counts, fill=geneBiotype)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = getPalette(colourCount))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Distribution of STAR counts per sample")
# ggsave("Results/Results1/DistrOfSTARcounts.png",p,width = 8,height = 7)


##HISAT2
##Note: uncomment ggsave() if you want to save plot

ggplot(data=dfHISAT2, aes(x=Sample, y=Percent_of_Counts, fill=geneBiotype)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values = getPalette(colourCount))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Distribution of HISAT2 counts per sample")
# ggsave("Output.png",p,width = 8,height = 7)


###################################################
## The number of identified protein-coding genes ##
###################################################

nidentifiedSTAR <- data.frame(countDataSTAR[1,], check.names = F)
nidentifiedHISAT2 <- data.frame(countDataHISAT2[1,], check.names = F)

for (i in 1:length(biotypes)) {
  n1<-which(countDataSTAR$gene_biotype==biotypes[i])
  n2<-which(countDataHISAT2$gene_biotype==biotypes[i])
  
  c1 <- biotypes[i]
  c2 <- biotypes[i]
  
  for (j in 2:ncol(countDataHISAT2)) {
    c1<-c(c1,length(which(countDataSTAR[n1,j]>0)))
    c2<-c(c2,length(which(countDataHISAT2[n2,j]>0)))

  }

  names(c1)<-colnames(countDataHISAT2)
  names(c2)<-colnames(countDataSTAR)
  
  nidentifiedSTAR<-rbind(nidentifiedSTAR,c1)
  nidentifiedHISAT2<-rbind(nidentifiedHISAT2,c2)
}

nidentifiedSTAR<-nidentifiedSTAR[-1,]
nidentifiedHISAT2<-nidentifiedHISAT2[-1,]

Tool<-rep("HISAT2",ncol(nidentifiedHISAT2)-1)
Samples<-colnames(nidentifiedHISAT2)[-1]
Number_of_protein_coding_genes<-data.matrix(nidentifiedHISAT2[which(nidentifiedHISAT2$gene_biotype=="protein_coding"),])[-1]

df1<-data.frame(Tool=Tool,Samples=Samples,Number_of_protein_coding_genes=Number_of_protein_coding_genes)


Tool<-rep("STAR",ncol(nidentifiedSTAR)-1)
Samples<-colnames(nidentifiedSTAR)[-1]
Number_of_protein_coding_genes<-data.matrix(nidentifiedSTAR[which(nidentifiedSTAR$gene_biotype=="protein_coding"),])[-1]

df2<-data.frame(Tool=Tool,Samples=Samples,Number_of_protein_coding_genes=Number_of_protein_coding_genes)

df<-rbind(df1,df2)

##Note: uncomment ggsave() function if you want to save the plot
ggplot(df, aes(factor(Samples), Number_of_protein_coding_genes, fill = Tool)) + 
  geom_bar(width=0.6,stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Number of detected protein coding genes")
# ggsave("Output.png",p,width = 9,height = 7)
