require(ggplot2)


###################################
## Library Diversity/Fisher test ##
###################################




############################################
## 1. MAPPING: Unique, Multiple, Unmapped ##
############################################

INPUTDIR <- "saminfo/"
files <- list.files(path = INPUTDIR, pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
files <- files[-which(files=="metaInfoMapping.txt")]
metaInfoMapping <- read.table(paste(INPUTDIR,"metaInfo.txt",sep=""), header = T,
                       check.names = F)

dat <- data.frame(Sample=NA, Percent=NA, Number=NA, Alignment=NA)
for (i in 1:length(files)) {
  samoutput <- read.table(paste(INPUTDIR,files[i],sep=""), 
                          header = FALSE,
                          sep="\t", fill=T)
  
  n <- as.numeric(as.character(samoutput[5,2]))
  
  nunique <- as.numeric(as.character(samoutput[8,2]))
  percent_unique <- nunique/n
  
  nmulti <- as.numeric(as.character(samoutput[23,2])) + 
            as.numeric(as.character(samoutput[25,2]))
  
  percent_multi <- nmulti/n
  
  nunmapped <-  as.numeric(as.character(samoutput[28,2]))+
                as.numeric(as.character(samoutput[30,2]))+
                as.numeric(as.character(samoutput[32,2]))
  
  percent_unmapped <- nunmapped/n
  
  filename <- unlist(strsplit(files[i], split=".", fixed=T))[1]
  
  newdata <- data.frame(Sample=rep(filename,3), 
                        Percent=c(percent_unique,percent_multi,percent_unmapped),
                        Number=c(nunique, nmulti,nunmapped),
                        Alignment=c("Unique alignment","Multiple alignment","Unmapped"))
  
  dat <- rbind(dat, newdata)
}

dat <- dat[-1,]

##Plot to see alignment for paired FFPE-FF samples
ggplot(data=dat, aes(x=Sample, y=Percent, fill=Alignment)) +
  geom_bar(stat="identity", position=position_dodge(), width = 0.7) + 
  theme(axis.text.x=element_text(angle=90, hjust=1))

# ggsave("../saminfo/Alignment.png", width = 9, height = 5)


##Fisher test
fisher_mapping <- data.frame(Sample=NA, 
                  Total=NA,
                  Uniquely_mapped=NA,
                  Not_uniquely_mapped_or_unmapped=NA,
                  Fisher_Exact_test_unique=NA,
                  Multiple_mapped	=NA,
                  Not_multiple_mapped_or_unmapped = NA,
                  Fisher_Exact_test_multiple =NA,
                  Unmapped=NA,
                  Mapped_unique_and_multiple = NA,
                  Fisher_Exact_test_unmapped = NA)



for (i in nrow(metaInfoMapping)) {

  samoutput_FF <- read.table(paste(INPUTDIR,metaInfoMapping$FF_filename[i],sep=""), 
                          header = FALSE,
                          sep="\t", fill=T)
  samoutput_FFPE <- read.table(paste(INPUTDIR,metaInfoMapping$FFPE_filename[i],sep=""), 
                             header = FALSE,
                             sep="\t", fill=T)
  
  n_FF <- as.numeric(as.character(samoutput_FF[5,2]))
  n_FFPE <- as.numeric(as.character(samoutput_FFPE[5,2]))
  
  nunique_FF <- as.numeric(as.character(samoutput_FF[8,2]))
  nunique_FFPE <- as.numeric(as.character(samoutput_FFPE[8,2]))
  
  nmulti_FF <- as.numeric(as.character(samoutput_FF[23,2])) + 
    as.numeric(as.character(samoutput_FF[25,2]))
  nmulti_FFPE <- as.numeric(as.character(samoutput_FFPE[23,2])) + 
    as.numeric(as.character(samoutput_FFPE[25,2]))
  
  nunmapped_FF <-  as.numeric(as.character(samoutput_FF[28,2]))+
    as.numeric(as.character(samoutput_FF[30,2]))+
    as.numeric(as.character(samoutput_FF[32,2]))

  nunmapped_FFPE <-  as.numeric(as.character(samoutput_FFPE[28,2]))+
    as.numeric(as.character(samoutput_FFPE[30,2]))+
    as.numeric(as.character(samoutput_FFPE[32,2]))
  
  mtrx_unique <- matrix(c(nunique_FF, n_FF-nunique_FF, 
                          nunique_FFPE, n_FFPE-nunique_FFPE), 
                 nrow=2,
                 ncol=2,
                 byrow = T)
  mtrx_multiple <- matrix(c(nmulti_FF, n_FF-nmulti_FF, 
                            nmulti_FFPE, n_FFPE-nmulti_FFPE), 
                        nrow=2,
                        ncol=2,
                        byrow = T)
  mtrx_unmapped <- matrix(c(nunmapped_FF, n_FF-nunmapped_FF, 
                            nunmapped_FFPE, n_FFPE-nunmapped_FFPE), 
                          nrow=2,
                          ncol=2,
                          byrow = T)
  
  res1 <- fisher.test(mtrx_unique, alternative = "greater")
  res2 <- fisher.test(mtrx_multiple, alternative = "less")
  res3 <- fisher.test(mtrx_unmapped, alternative = "less")
  
  filename_FF <- unlist(strsplit(as.character(metaInfoMapping$FF_filename[i]), split=".", fixed=T))[1]
  filename_FFPE <- unlist(strsplit(as.character(metaInfoMapping$FFPE_filename[i]), split=".", fixed=T))[1]
  
  newdata <- data.frame(Sample=c(filename_FF, filename_FFPE), 
                        Total = c(n_FF,n_FFPE),
                        Uniquely_mapped=c(nunique_FF,nunique_FFPE),
                        Not_uniquely_mapped_or_unmapped=c(n_FF-nunique_FF,n_FFPE-nunique_FFPE),
                        Fisher_Exact_test_unique=c(res1$p.value,res1$p.value),
                        Multiple_mapped	=c(nmulti_FF,nmulti_FFPE),
                        Not_multiple_mapped_or_unmapped = c(n_FF-nmulti_FF,n_FFPE-nmulti_FFPE),
                        Fisher_Exact_test_multiple =c(res2$p.value,res2$p.value),
                        Unmapped=c(nunmapped_FF,nunmapped_FFPE),
                        Mapped_unique_and_multiple = c(n_FF-nunmapped_FF,n_FFPE-nunmapped_FFPE),
                        Fisher_Exact_test_unmapped = c(res3$p.value,res3$p.value))
  
  fisher_mapping <- rbind(fisher_mapping, newdata)
}

fisher_mapping <- fisher_mapping[-1,] #result
# write.csv(fisher_mapping, "Mapping_FisherTest.csv")


###################
## 2. DUPLICATES ##
###################

duplicates_dat <- read.table("duplicates/Duplicates.txt",
                             header=T)

metaInfoDup <- read.table("duplicates/MetaInfoDup.txt", header = T)

fisher_dup <- data.frame(Sample=NA, 
                  Total=NA,
                  Duplicates=NA,
                  Unique=NA,
                  Fisher_Exact_test=NA)

for (i in 1:nrow(metaInfoDup)) {
  
  FFsample <- as.character(metaInfoDup$FF[i])
  FFPEsample <- as.character(metaInfoDup$FFPE[i])
  
  mtrx_dup <- matrix(c(duplicates_dat$Duplicate_Reads[which(duplicates_dat$Category==FFsample)], duplicates_dat$Unique_Reads[which(duplicates_dat$Category==FFsample)],
                       duplicates_dat$Duplicate_Reads[which(duplicates_dat$Category==FFPEsample)], duplicates_dat$Unique_Reads[which(duplicates_dat$Category==FFPEsample)]), 
                        nrow=2,
                        ncol=2,
                        byrow = T)
  
  n_FF <- duplicates_dat$Unique_Reads[which(duplicates_dat$Category==FFsample)]+duplicates_dat$Duplicate_Reads[which(duplicates_dat$Category==FFsample)]
  n_FFPE <- duplicates_dat$Unique_Reads[which(duplicates_dat$Category==FFPEsample)]+duplicates_dat$Duplicate_Reads[which(duplicates_dat$Category==FFPEsample)]
  
  res <- fisher.test(mtrx_dup, alternative = "less")
  
  newdata <- data.frame(Sample=c(FFsample, FFPEsample), 
                        Total = c(n_FF,n_FFPE),
                        Duplicates=c(duplicates_dat$Duplicate_Reads[which(duplicates_dat$Category==FFsample)],duplicates_dat$Duplicate_Reads[which(duplicates_dat$Category==FFPEsample)]),
                        Unique=c(duplicates_dat$Unique_Reads[which(duplicates_dat$Category==FFsample)],duplicates_dat$Unique_Reads[which(duplicates_dat$Category==FFPEsample)]),
                        Fisher_Exact_test=c(res$p.value,res$p.value))
  
  fisher_dup <- rbind(fisher_dup, newdata)
}

fisher_dup <- fisher_dup[-1,] #result
# write.csv(fisher_dup, "Duplicates_Fishertest.csv")



##########################
## 3. Read Distribution ##
##########################

##Fisher test

dat_readDist <- read.table("read_distribution/rseqc_read_distribution_paired.txt", header=T,
                           check.names = F)

metaInfo_readDist <- read.table("read_distribution/MetaInfo_ReadDist.txt", header = T)

fisher_readDistr <- data.frame(Sample=NA,
                        Total=NA,
                        Exonic=NA,
                        Non_exonic=NA,
                        Fisher_Exact_test_exonic=NA,
                        Intronic=NA,
                        Non_intronic=NA,
                        Fisher_Exact_test_intronic=NA,
                        Intergenic=NA,
                        Non_intergenic=NA,
                        Fisher_Exact_test_intergenic=NA
                        )

for (i in 1:nrow(metaInfo_readDist)) {

  sample_FF <- as.character(metaInfo_readDist$FF[i])
  sample_FFPE <- as.character(metaInfo_readDist$FFPE[i])
  
  i_FF <- which(dat_readDist$Category==sample_FF)
  i_FFPE <- which(dat_readDist$Category==sample_FFPE)
  
  intronic_FF <- dat_readDist$Introns[i_FF]
  intronic_FFPE <- dat_readDist$Introns[i_FFPE]
  
  exonic_FF <- dat_readDist$CDS_Exons[i_FF]+dat_readDist$`5'UTR_Exons`[i_FF]+dat_readDist$`3'UTR_Exons`[i_FF]
  exonic_FFPE <- dat_readDist$CDS_Exons[i_FFPE]+dat_readDist$`5'UTR_Exons`[i_FFPE]+dat_readDist$`3'UTR_Exons`[i_FFPE]
  
  intergenic_FF <- dat_readDist$TSS_up_10kb[i_FF]+dat_readDist$TES_down_10kb[i_FF]+dat_readDist$Other_intergenic[i_FF]
  intergenic_FFPE <- dat_readDist$TSS_up_10kb[i_FFPE]+dat_readDist$TES_down_10kb[i_FFPE]+dat_readDist$Other_intergenic[i_FFPE]
  
  n_FF <- intronic_FF+exonic_FF+intergenic_FF
  n_FFPE <- intronic_FFPE+exonic_FFPE+intergenic_FFPE
  
  
  mtrx_exonic <- matrix(c(exonic_FF,exonic_FFPE,
                          n_FF-exonic_FF,n_FFPE-exonic_FFPE), 
                     nrow=2,
                     ncol=2,
                     byrow = T)
  
  res1 <- fisher.test(mtrx_exonic, alternative = "greater")
  
  mtrx_intronic <- matrix(c(intronic_FF,intronic_FFPE,
                          n_FF-intronic_FF,n_FFPE-intronic_FFPE), 
                        nrow=2,
                        ncol=2,
                        byrow = T)
  
  res2 <- fisher.test(mtrx_intronic, alternative = "less")
  
  mtrx_intergenic <- matrix(c(intergenic_FF,intergenic_FFPE,
                            n_FF-intergenic_FF,n_FFPE-intergenic_FFPE), 
                          nrow=2,
                          ncol=2,
                          byrow = T)
  
  res3 <- fisher.test(mtrx_intergenic, alternative = "less")
  
  newrow <- data.frame(Sample=c(sample_FF,sample_FFPE),
             Total=c(n_FF,n_FFPE),
             Exonic=c(exonic_FF,exonic_FFPE),
             Non_exonic=c(n_FF-exonic_FF,n_FFPE-exonic_FFPE),
             Fisher_Exact_test_exonic=c(res1$p.value,res1$p.value),
             Intronic=c(intronic_FF,intronic_FFPE),
             Non_intronic=c(n_FF-intronic_FF,n_FFPE-intronic_FFPE),
             Fisher_Exact_test_intronic=c(res2$p.value,res2$p.value),
             Intergenic=c(intergenic_FF,intergenic_FFPE),
             Non_intergenic=c(n_FF-intergenic_FF,n_FFPE-intergenic_FFPE),
             Fisher_Exact_test_intergenic=c(res3$p.value,res3$p.value))
             
  fisher_readDistr <- rbind(fisher_readDistr,newrow)

}
fisher_readDistr <- fisher_readDistr[-1,] #result
# write.csv(readDistr,"../read_distribution/readDistr.csv")


##PLOT

readDistrPlot <- data.frame(Sample=NA,
                            Percent=NA,
                            Region=NA,
                            Pair=NA)

for (i in 1:nrow(metaInfo_readDist)) {
  
  sample_FF <- as.character(metaInfo_readDist$FF[i])
  sample_FFPE <- as.character(metaInfo_readDist$FFPE[i])
  
  i_FF <- which(dat_readDist$Category==sample_FF)
  i_FFPE <- which(dat_readDist$Category==sample_FFPE)
  
  intronic_FF <- dat_readDist$Introns[i_FF]
  intronic_FFPE <- dat_readDist$Introns[i_FFPE]
  
  exonic_FF <- dat_readDist$CDS_Exons[i_FF]+dat_readDist$`5'UTR_Exons`[i_FF]+dat_readDist$`3'UTR_Exons`[i_FF]
  exonic_FFPE <- dat_readDist$CDS_Exons[i_FFPE]+dat_readDist$`5'UTR_Exons`[i_FFPE]+dat_readDist$`3'UTR_Exons`[i_FFPE]
  
  intergenic_FF <- dat_readDist$TSS_up_10kb[i_FF]+dat_readDist$TES_down_10kb[i_FF]+dat_readDist$Other_intergenic[i_FF]
  intergenic_FFPE <- dat_readDist$TSS_up_10kb[i_FFPE]+dat_readDist$TES_down_10kb[i_FFPE]+dat_readDist$Other_intergenic[i_FFPE]
  
  n_FF <- intronic_FF+exonic_FF+intergenic_FF
  n_FFPE <- intronic_FFPE+exonic_FFPE+intergenic_FFPE

  newrow_FF <- data.frame(Sample=rep(sample_FF,3),
                       Percent=c(round((exonic_FF/n_FF)*100),round((intronic_FF/n_FF)*100),round((intergenic_FF/n_FF)*100)),
                       Region=c("Exonic","Intronic","Intergenic"),
                       Pair = rep(paste("pair",i),3)
                       )
  
  newrow_FFPE <- data.frame(Sample=rep(sample_FFPE,3),
                        Percent=c(round((exonic_FFPE/n_FFPE)*100),round((intronic_FFPE/n_FFPE)*100),round((intergenic_FFPE/n_FFPE)*100)),
                        Region=c("Exonic","Intronic","Intergenic"),
                        Pair = rep(paste("pair",i),3))
  
  
  readDistrPlot  <- rbind(newrow_FF ,newrow_FFPE)
  
  ggplot(data=readDistrPlot, aes(x=Sample, y=Percent, fill=Region)) +
    geom_bar(stat="identity", position=position_dodge(), width = 0.5) + facet_wrap(~Pair)+
    scale_fill_manual(values=c("#CC6666", "#9999CC", "#66CC99"))
  # ggsave("plot.png", width = 6, height = 3)
}