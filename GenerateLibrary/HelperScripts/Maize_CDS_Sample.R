library(tidyverse)
library(readr)
library(readxl)
library(dplyr)

maize_CDS <- read.table("RR_maize_V4_CDS_170_length.fa",header=FALSE)
library(stringr)
# filter out sequences with a BsaI or Esp3I site (the cloning adapter adds a G in front of the sequence)
RE.sites <- c('(G|C|^)GTCTC', '(G|^)AGAC(C|G)')
tmpy <- maize_CDS$V2[! grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), maize_CDS$V2)]
maize_CDS_noRE <- maize_CDS %>% filter(V2 %in% tmpy) %>% mutate(GC=((str_count(V2,"G"))+(str_count(V2,"C")))/str_length(V2)) %>% filter(GC < .55)

set.seed(10000)
maize_CDS_500_sample <- sample_n(maize_CDS_noRE,589)
maize_CDS_500_sample %>% mutate(name=paste(V1,"CDS",sep="_")) %>% select(name,V2) %>% write_delim('Maize_CDS_589_sample.fa', col_names = FALSE, delim = "\n")
