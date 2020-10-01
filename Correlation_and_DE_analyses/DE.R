##DESCRIPTION:
##Differential expression analysis for both FFPE and fresh-frozen samples
##The following groups are compared: 
##FFPE_od versus FFPE_control; and FF_od versus FF_control
##DESeq2 package is used for the DE measurements

##Include packages (install if necessary)
list.of.packages <- c("DESeq2","RColorBrewer","gplots","ggplot2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library("DESeq2")
library("RColorBrewer")
library("gplots")
library("ggplot2")

##define custom functions:
##DEanalysis (cond1, cond2)
##PCAplot (vsd)
##NumberOfProtenCodingGenes (res_FF,res_FFPE)
##Heatmap_TAAs (vsd)
source("DE_functions.R") 

source("biomart.R") ##Load biomart db containing information about canine gene ids,gene names, gene biotypes etc.. The seult is wtitten into the "BM" data.frame


################
## DATA INPUT ##
################

countData<-read.table('INPUT/InputCountsFC_STAR.txt', header=TRUE, check.names=FALSE)

##remove samples which failed QC in our RNA-seq (see text of the thesis)
id_rm<-c()
id_rm<-c(id_rm,which(colnames(countData)=="50315_FFPE_control"))
id_rm<-c(id_rm,which(colnames(countData)=="60462_FFPE_control"))
id_rm<-c(id_rm,which(colnames(countData)=="44331_FFPE_od"))
countData<-countData[,-id_rm]


rownames(countData)<-countData$GeneID
countData<-countData[, -c(1:2)] ##remove meatadata not needed for the DE analysis

##Get sample groups (FFPE_od/FFPE_control/FF_od/FF_control)
##Each sample name has the following structure: {patient id}_[FFPE/FF]_[od/control]
##The part "[FFPE/FF]_[od/contro]" implies the sample group
groups <- c()
for (i in c(1:length(colnames(countData)))) {
  conservType<-unlist(strsplit(colnames(countData)[i],"_",fixed=T))[2]
  condition<-unlist(strsplit(colnames(countData)[i],"_",fixed=T))[3]
  group<-paste(conservType,condition,sep = "_")
  groups <- c(groups,group)
}

colData<-data.frame(condition=factor(groups)) ##samples information


####################################
## 1. DESeq2 EXPRESSSION ANALYSIS ##
####################################


dds<-DESeqDataSetFromMatrix(countData, colData, formula(~condition))
dds<-DESeq(dds, betaPrior = TRUE) ##betaPrior = TRUE makes lfcSrhrinking
#res <-results(dds)
vsd <-vst(dds, blind=TRUE) #Variance Stabilizing Transformation

##Note: uncomment png() and dev.off() functions in case you want to save the plot

# png("PCAnew.png", width = 800, height = 800)
pca <- PCAplot(vsd) ##PCA plot
plot(pca)
# dev.off()

##To see what are the pairs P1-P8 --> check pca[["data"]]

################################################
## 2. DESeq2 DIFFERENTIAL EXPRESSION ANALYSIS ##
################################################

counts <- counts(dds, normalized=T) ##Normalized DESeq2 counts
metaData<-as.data.frame(dds@colData) ##metadata

res_FF <- DEanalysis("FF_od","FF_control") #DE result for fresh-frozen
res_FFPE <- DEanalysis("FFPE_od","FFPE_control") #DE result for FFPE

#MAplot
plotMA(results(dds, contrast=c("condition","FF_od","FF_control")), alpha=0.01, main="MAplot fresh-frozen")
plotMA(results(dds, contrast=c("condition","FFPE_od","FFPE_control")), alpha=0.01, main="MAplot FFPE")

#Volcano plot

##NOTE:Uncomment png() and dev.off() functions for saving the plot in function VolcanoPlot()

VolcanoPlot(res_FF, "fresh-frozen")
VolcanoPlot(res_FFPE, "FFPE")

####################################
## 3. Filtering of the DE results ##
####################################


#Setting the threshold of significance. Leave only DE genes with FDR<0.01, |logFC|>1 and baseMean>10
FDR_lim=0.01
logFC_lim=1
baseMean_lim=10

#fresh-frozen
res_FF<-res_FF[which(res_FF$padj<FDR_lim &
                     abs(res_FF$log2FoldChange)>logFC_lim &
                     res_FF$baseMean>baseMean_lim),]


#FFPE
res_FFPE<-res_FFPE[which(res_FFPE$padj<FDR_lim &
                         abs(res_FFPE$log2FoldChange)>logFC_lim &
                         res_FFPE$baseMean >baseMean_lim),]


##**UNCOMMENT TO WRITE THE RESULTS INTO THE FILES**
# write.csv(res_FFPE,
#           file="Results_FFPE_od_versus_FFPE_control.csv",row.names = F)
# write.csv(res_FF,
#           file="Results_FF_od_versus_FF_control.csv",row.names = F)



##Build a bar plot with number of DE protein-coding genes. Comparison of results of DE FFPE and DE FF analyses 
NumberOfProtenCodingGenes (res_FF,res_FFPE)

##Building of heatmap with all defined TAA candidates (see "TAAs.txt")  
Heatmap_TAAs(vsd)
