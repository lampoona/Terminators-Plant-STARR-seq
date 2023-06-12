library(tidyverse)
library(readr)
library(readxl)
library(dplyr)

Genes <- read.table("maize_V4_CDS_clean_sorted.bed", header=FALSE)
Genes$length <- Genes$V3-Genes$V2

#take the genes greather than 170 in length and just take out the first 170
Genes170 <- Genes %>% filter(length >= 170) 
Genes170$middle <- floor((Genes170$V3+Genes170$V2)/2) 
Genes170$V2 <- Genes170$middle - 85
Genes170$V3 <- Genes170$middle + 85
Genes170 <- Genes170 %>% select(-length,-middle) 
Genes170_1 <- Genes170 %>% arrange(by=V1,V2)
write.table(Genes170_1, file="maize_V4_CDS_170_length.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
