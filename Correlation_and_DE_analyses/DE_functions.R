##FUNCTION NEEDED FOR THE DE ANALYSIS

#################################
## DESEQ2 DE ANALYSIS FUNCTION ##
#################################

DEanalysis <- function(cond1, cond2) { ##DE analysis: "cond1" versus "cond2"
  
  ##DE analysis between cond1 and cond2 conditions
  counts_od<-counts[,rownames(metaData)[which(metaData$condition==cond1)]]
  counts_control<-counts[,rownames(metaData)[which(metaData$condition==cond2)]]
  res<-as.data.frame(results(dds, contrast=c("condition",cond1,cond2))) ##DE analysis
  res <- cbind(res,counts_od,counts_control)
  res <- data.frame(GeneID = row.names(res), res, check.names = FALSE) 
  
  ##Add to the DE results data frame information about gene name, gene biotype and description of gene
  res[,"GeneName"]<-NA #Initialization of new column
  res[,"GeneBiotype"]<-NA #Initialization of new column
  res[,"GeneDescription"]<-NA #Initialization of new column
  for (i in 1:length(res$GeneID)) {
    geneid<-res$GeneID[i]
    res$GeneName[i]<-BM$external_gene_name[which(BM$ensembl_gene_id==geneid)]
    res$GeneBiotype[i]<-BM$gene_biotype[which(BM$ensembl_gene_id==geneid)]
    res$GeneDescription[i]<-genedescription<-BM$description[which(BM$ensembl_gene_id==geneid)]
  }
  
  
  res<-res[,c(1,ncol(res)-2,ncol(res)-1,ncol(res),2:(ncol(res)-3))] ##re-order columns
  
  ##Re-calculating baseMean (by default, calculated for all 4 input groups)
  for (i in 1:nrow(res)) {
    base_mean <- mean(as.numeric(res[i,12:ncol(res)]))
    res$baseMean[i]<- base_mean
  }
  
  return(res) #DE result
}




#############################################################
## BUILDING OF PCA PLOT WITH MARKED PAIRED FFPE-FF SAMPLES ##
#############################################################

PCAplot <- function(vsd) {

  pairedSamples <- read.table("INPUT/Paired_FFPEvsFF.txt", header = T, check.names = F)
  
  pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  pcaData["paired"] <- NA
  
  #Marking matching FFPE-FF samples. Unpaired samples remain NA
  for (i in 1:length(pcaData$name)) {
    for (j in 1:nrow(pairedSamples)) {
      if (pcaData$name[i] == pairedSamples$FFPE[j]| pcaData$name[i] == pairedSamples$FF[j]) pcaData$paired[i] <- paste("P",j,sep = "") 
    }
  }
  
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  ##The plot can be saved in the PNG format**
  ##Change "OUTPUT" to desired name**
  ##Then, uncomment png() and dev.off() fucntions**
  
  p <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed() + 
    geom_line(data = pcaData[!is.na(pcaData$paired), ], 
              aes(group = paired),
              linetype = "dashed",
              color="grey", size=1) + 
    geom_text(aes(label=pcaData$paired,vjust=2))  
  return (p)

}



#####################################################################
## FUNCTION BUILDS BAR PLOT WITH NUMBER OF DE PROTEIN-CODING GENES ##
## DETECTED IN BOTH FFPE DE AND FF DE ANALYSES FOR COMPARISON      ##
#####################################################################

NumberOfProtenCodingGenes <- function (res_FF,res_FFPE) {
  
  res_FF_protein_coding <- res_FF[which(res_FF$GeneBiotype=="protein_coding"),] ##protein-coding genes only
  res_FFPE_protein_coding <- res_FFPE[which(res_FFPE$GeneBiotype=="protein_coding"),] ##protein-coding genes only 
  
  #The number of up- and down- regulated protein-coding genes in FF samples
  n_FF_up <- length(which(res_FF_protein_coding$log2FoldChange>0)) #up-regulated
  n_FF_down <- length(which(res_FF_protein_coding$log2FoldChange<0)) #down-regulated
  
  #The number of up- and down- regulated protein-coding genes in FFPE samples
  n_FFPE_up <- length(which(res_FFPE_protein_coding$log2FoldChange>0)) #up-regulated
  n_FFPE_down <- length(which(res_FFPE_protein_coding$log2FoldChange<0)) #down-regulated
  
  DEdat <- data.frame(Group=c("fresh-frozen","fresh-frozen","FFPE","FFPE"),
                      Regulation_in_od=c("down","up","down","up"),
                      Number_of_protein_genes=c(n_FF_down,n_FF_up,n_FFPE_down,n_FFPE_up))
  
  #build the plot
  ggplot(data=DEdat, aes(x=Group, y=Number_of_protein_genes, fill=Regulation_in_od)) +
    geom_bar(stat="identity", position=position_dodge(), width = 0.5)+
    scale_fill_manual(values = c("#0392cf","#ee4035")) 
}


VolcanoPlot <- function(res, label) {
  
  ##NOTE:Uncomment png() and dev.off() functions for saving the plot**
  
  # png(filename = "Results/DEgenes/VolcanoPlotFF.png", width = 1500, height = 1500)
  par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
  
  with(res, plot(log2FoldChange, -log10(padj), pch=20, main=paste("Volcano plot",label, sep=" "), cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
  with(subset(res, padj<0.01 & abs(log2FoldChange)>1 & baseMean>10), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
  with(subset(res, -log10(padj)>6 & abs(log2FoldChange)>1), text(log2FoldChange, -log10(padj), labels=subset(res$GeneName,-log10(res$padj)>6 & abs(res$log2FoldChange)>1), cex=0.5, pos=2))
  
  abline(v=0, col="black", lty=3, lwd=1.0)
  abline(v=-1, col="black", lty=4, lwd=2.0)
  abline(v=1, col="black", lty=4, lwd=2.0)
  abline(h=-log10(max(res$padj[res$padj<0.01], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
  # dev.off()
}


#############################################################
## BUILDING OF HEATMAP WITH ALL DEFINED TAA CANDIDATES     ##
## NOTE: the TAA candidates have been defined using also   ##
## some third-party applications e.g. IPA.                 ##
## The list of the TAA candidates is defined in "TAAs.txt" ##
#############################################################

Heatmap_TAAs <- function(vsd) {
  
  vsd <- assay(vsd)
  
  TAAcandidates <- read.table("INPUT/TAAs.txt", header = T,
                              stringsAsFactors = F)
  TAAcandidates <- TAAcandidates$TAA_geneName
  
  
  TAAinfo <- data.frame(Genename="",
                        Geneid="")
  for (i in 1:length(TAAcandidates)) {
    genename <- TAAcandidates[i]
    geneid <- BM$ensembl_gene_id[which(BM$external_gene_name==genename)]
    newrow <- data.frame(Genename=genename,
                         Geneid=geneid)
    TAAinfo <- rbind(TAAinfo,newrow)
    
  }
  TAAinfo <- TAAinfo[-1,]
  
  select <- c()
  for(i in 1:length(TAAinfo$Geneid)) {
    geneid <- as.character(TAAinfo$Geneid[i])
    find <- which(rownames(vsd)==geneid)
    select <- c(select,find)
    rownames(vsd)[find] <- as.character(TAAinfo$Genename[i])
  }
  
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  
  ##NOTE:Uncomment pdf() and dev.off() functions for saving the heatmap**
  
  #pdf(paste("Results/heatmap/","TAAcandidates.pdf",sep=""),width = 10,height = 13)
  heatmap.2(vsd[select,], col = hmcol, trace="none", margin=c(10,6),
            labCol=colnames(vsd), cexRow = 1/log10(length(select)))
  #dev.off()
}
