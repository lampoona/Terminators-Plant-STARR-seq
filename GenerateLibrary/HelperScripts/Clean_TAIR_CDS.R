library(tidyverse)
library(readr)
library(readxl)
library(dplyr)

Genes <- read.table("/net/gs/vol1/home/kbubb/queitschlab/TAIR_RELEASE_GENOMES/TAIR10/TAIR10_GFF3_CDS.gff", header=FALSE)
Genes$length <- Genes$V5-Genes$V4
Genes170 <- Genes %>% filter(length >= 170)
Genes170_2 <- Genes170 %>% separate(V9,into = c("Parent","name"), sep = "=") %>% separate(name,c("name", "blah"), sep=",") %>% select(-Parent,-blah,-V2,-V8,-V3,length) #%>% relocate(name, .after=V5)
data <- Genes170_2[,c(1,2,3,6,4,5)]
colnames(data) <- c("V1","V2","V3","V4","V5","V6")

data$middle <- round((data$V3+data$V2)/2) 
data$V2 <- data$middle - 85
data$V3 <- data$middle + 85
data <- data %>% select(-middle)
foo <- data %>% filter(V1 != "ChrM") %>% filter(V1 != "ChrC")
write.table(foo, file="TAIR10_CDS_170_plus.bed", col.names = FALSE, row.names = FALSE, quote = FALSE, sep="\t")
