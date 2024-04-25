#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

library(dplyr)
library(tidyr)


peaker = function(x){
  flank=3  # how many to look at on either side
  peaks=rep(0,length(x))
  vector=rep(0,length(x))
  total=sum(x)
  for(i in c((flank+1):(length(x)-flank))){
    if(x[i] == max(x[c((i-flank):(i+flank))])){ # if pos is higher than the 4 pos on either side
      if(x[i] > (0.01 * total)){ # if pos accounts for more than 1% of total positions
        peaks[i]=1
        vector[i]=1
      }
    }
  }
  for(i in c(219:222)){
         if(x[i]> (0.01 * total)){ 
                peaks[222]=1 
                vector[222]=1
                }
        }
        df <- data.frame("peaks"=peaks,"place"=vector)
  return(df)
}

blah <- read.table(args[1],header=FALSE)
colnames(blah) <- c("RN","UMI","gene","seq")
gene = blah %>% pull(gene) %>% unique() %>% as.character()
blah <- blah %>% mutate(length=nchar(as.character(seq)))# %>% filter(length >= 10) %>% select(UMI,length) 
dedup <- blah %>% count(UMI,length) %>% select(-n)
trim <- dedup %>% filter(length > 10) 

if(length(unique(trim$length))==1) {
  save <- trim %>% count(length) %>% as.data.frame()
  colnames(save) <- c("peak","counts")
  write.table(save, paste(gene,"/",gene,"_PACs.tsv",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 
}

#need to add a read count limit for peak calling. Perhaps less than 10 reads is not good?

if(length(unique(trim$length))>1) {
  same <- max(trim$length)-min(trim$length)
  counts <- hist(trim$length,breaks=same)$counts
  breaks <- hist(trim$length, breaks=same)$breaks 
  names(counts) <- breaks[-length(breaks)]


  df <- as.data.frame(counts)
  df <- tibble::rownames_to_column(df, "position")
  df <- df %>% mutate(position = as.numeric(position) + 1)


  dummy <- data.frame(position=seq(1,222,1))
  final <- dummy %>% left_join(df,by='position')
  final[is.na(final)] = 0


#write.table(final,paste(gene,"/",gene,"_hist_PACs.tsv",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 
	
  pdf(paste(gene,"/",gene,".pdf", sep=""))
  plot(final$position, final$counts, type="h",ylab="counts",xlab="cleavage site",main=gene)
  peaks1 = peaker(final$counts)
  points(final$position[which(peaks1$place == 1)],final$counts[which(peaks1$place == 1)], col="red")
  dev.off()

  output <- data.frame(peak=final$position[which(peaks1$place == 1)],counts=final$counts[which(peaks1$peaks == 1)])
	
  output2 <- data.frame(position=final$position[which(peaks1$place == 1)],counts=final$counts[which(peaks1$peaks == 1)])
  output2 <- output2 %>% mutate(peaker=counts)
  final <- final %>% full_join(output2,by=c("position","counts")) %>% replace(is.na(.), 0)
  write.table(final,paste(gene,"/",gene,"_hist_PACs.tsv",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  
  write.table(output, paste(gene,"/",gene,"_PACs.tsv",sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE) 
}
