require("ggplot2")


#####################################################
## Plot innre-distance of matching FFPE-FF samples ##
#####################################################

INPUTDIR <- "inner_distance/"

files <- list.files(path = INPUTDIR, pattern = NULL, all.files = FALSE,
                    full.names = FALSE, recursive = FALSE,
                    ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

files <- files[-which(files=="MetaInfo.txt")]

meta <- read.table(paste(INPUTDIR,"MetaInfo.txt",sep=""),header = T)

for (i in 1:nrow(meta)) {
  
  freqFF <- read.table(paste(INPUTDIR,meta$FF[i],sep = ""))
  freqFFPE <- read.table(paste(INPUTDIR,meta$FFPE[i],sep = ""))
  
  x <- seq(-248,247, by=5)
  
  filename_FF <- unlist(strsplit(as.character(meta$FF[i]),split=".", fixed=TRUE))[1]
  filename_FFPE <- unlist(strsplit(as.character(meta$FFPE[i]),split=".", fixed=TRUE))[1]
  
  dat_FF <- data.frame(inner_distance=x, sample=filename_FF, frequency=freqFF$V3)
  dat_FFPE <- data.frame(inner_distance=x, sample=filename_FFPE, frequency=freqFFPE$V3)
  
  dat_FF <- dat_FF[which(dat_FF$frequency>0),]
  dat_FFPE <- dat_FFPE[which(dat_FFPE$frequency>0),]
  
  fragsize_FF=rep(dat_FF$inner_distance,times=dat_FF$frequency)
  frag_sd_FF = sd(fragsize_FF)
  frag_mean_FF = mean(fragsize_FF)
  frag_median_FF = median(fragsize_FF)
  p_FF <- hist(fragsize_FF,probability=T,breaks=100,freq=F,xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean_FF,";","SD=",frag_sd_FF),collapse=""),border="blue")

  fragsize_FFPE=rep(dat_FF$inner_distance,times=dat_FFPE$frequency)
  frag_sd_FFPE = sd(fragsize_FFPE)
  frag_mean_FFPE = mean(fragsize_FFPE)
  frag_median_FFPE = median(fragsize_FFPE)
  p_FFPE <- hist(fragsize_FFPE,probability=T,breaks=100,freq = F, xlab="mRNA insert size (bp)",main=paste(c("Mean=",frag_mean_FFPE,";","SD=",frag_sd_FFPE),collapse=""),border="blue")
  
  minx <- min(min(dat_FF$inner_distance),min(dat_FFPE$inner_distance)) 
  maxx <- max(max(dat_FF$inner_distance),max(dat_FFPE$inner_distance)) 
  maxy <- max(max(p_FF$density),max(p_FFPE$density))
  
  set.seed(42)
  
  # png(paste("../inner_distance/paired/plotsR/",unlist(strsplit(filename_FF, split = "_"))[1],".png", sep=""))
  plot( p_FF, col=rgb(0,0,1,1/4),freq = F, xlim=c(minx,maxx), ylim = c(0,maxy),
        main = paste(filename_FF, "vs.", filename_FFPE, sep=" "),
        xlab = "inner distance")  # first histogram
  plot( p_FFPE, col=rgb(1,0,0,1/4), freq=F, xlim=c(minx.maxx), ylim = c(0,maxy), add=T)  # second histogram
  # dev.off()

}
