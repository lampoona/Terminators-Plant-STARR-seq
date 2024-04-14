library(tidyverse)
library(readr)
library(readxl)
library(dplyr)

secondMax <-  function(x) {
  u <- unique(x)
  sort(u, decreasing = TRUE)[2L]
}
#File_S3_Maize_clean.txt  is a cleaned up version of supplemental table 3 in the  Jafar et al., 2019  tables (File S3) doi: 10.1534/g3.119.400196 
#https://gsajournals.figshare.com/articles/dataset/Supplemental_Material_for_Jafar_et_al_2019/8034695?file=14964725

data <- read.table("File_S3_Maize_clean.txt",header=TRUE)

### filter for protein-coding genes on the nuclear chromosomes and keep only the most abundant PAC for each gene ###
main_PACs <-  data %>%
  filter(Gene_biotype == 'protein_coding') %>%
  group_by(Chromosome, Strand, GeneID) %>%
  slice_max(
    order_by = PAC_reads,
    n = 1
  ) %>%
  mutate(
    helper = if_else(Strand == '+', PAC, -PAC)
  ) %>%
  slice_min(
    order_by = helper,
    n = 1
  ) %>%
  ungroup() %>%
  select(-helper) %>% mutate(place=rep("1"))


##get secondary PACs that have at least 30% of total PAC reads (of top two PACs)

secondary_PACs <- data %>% filter(Gene_biotype == 'protein_coding') %>% 
  group_by(TranscriptID) %>% 
  arrange(by=desc(PAC_reads),.by_group = TRUE) %>%
  slice_max(order_by = PAC_reads,n=2,with_ties = FALSE) %>%
  mutate(percent_of = PAC_reads/sum(PAC_reads)*100) %>%
  mutate(max = max(PAC_reads)) %>%
  mutate(helper = if_else(PAC_reads == secondMax(PAC_reads),TRUE,FALSE)) %>%
  filter(helper == TRUE) %>%
  filter(PAC_reads > 30) %>% 
  ungroup() %>% select(-max,-percent_of,-helper) %>% mutate(place = rep("2"))

Final_maize <- rbind(main_PACs,secondary_PACs) %>% group_by(GeneID) %>% ungroup()

#prepare final document as a bed file (170 total length)

Final_maize %>%
  mutate(
    Chromosome = if_else(Chromosome != "chr10",gsub('chr', 'chr0', Chromosome, fixed = TRUE),'chr10'),
    start = if_else(Strand == '+', PAC - 151, PAC - 21),
    end = if_else(Strand == '+', PAC + 19, PAC + 149),
    score = paste(PAC_reads, Cluster_reads, sep = '/')
  ) %>%
  select(Chromosome, start, end, TranscriptID, score, Strand,place) %>%
  arrange(Chromosome, start, TranscriptID) %>%
  mutate(TranscriptID = if_else(place == "2",paste(TranscriptID,"_2",sep = ''),TranscriptID)) %>% 
  select(-place) %>%
  write_tsv('Maize_Top_PACs_Final.bed', col_names = FALSE)







