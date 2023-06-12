library(tidyverse)
library(stringr)
library(dplyr)

TAIR_CDS <- read.table("RR_TAIR10_CDS_170_plus.fa",header=FALSE)
TAIR_CDS$V1<-substr(TAIR_CDS$V1,1,nchar(TAIR_CDS$V1)-2)

# filter out sequences with a BsaI or Esp3I site (the cloning adapter adds a G in front of the sequence)
RE.sites <- c('(G|C|^)GTCTC', '(G|^)AGAC(C|G)')
tmpy <- TAIR_CDS$V2[! grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), TAIR_CDS$V2)]
TAIR_CDS_noRE <- TAIR_CDS %>% filter(V2 %in% tmpy) %>% mutate(GC=((str_count(V2,"G"))+(str_count(V2,"C")))/str_length(V2)) %>% filter(GC < .45)
set.seed(10000)
TAIR_CDS_500_sample <- sample_n(TAIR_CDS_noRE,589)

TAIR_CDS_500_sample %>% mutate(name=paste(V1,"CDS",sep="_")) %>% select(name,V2) %>% write_delim('TAIR_CDS_589_sample.fa', col_names = FALSE, delim = "\n")
