library(tidyverse)
library(readr)
library(readxl)
library(dplyr)


secondMax <-  function(x) {
  u <- unique(x)
  sort(u, decreasing = TRUE)[2L]
}


arab_170 <- read.table("Arabidopsis_170.bed", header=FALSE)
arab_170$V4 <- substr(arab_170$V4,1,nchar(arab_170$V4)-2)

data_sense <- read_xlsx('Thomas2012_SuppTables.xlsx', sheet = 1, range = 'H14:J17987') %>%
  bind_rows(
    read_xlsx('Thomas2012_SuppTables.xlsx', sheet = 1, range = 'N14:P20606')
  )
data_antisense <- read_xlsx('Thomas2012_SuppTables.xlsx', sheet = 2, range = 'G14:J6491') %>%
  bind_rows(
    read_xlsx('Thomas2012_SuppTables.xlsx', sheet = 2, range = 'M14:P12839')
  ) %>%
  filter(as.integer(substr(type, nchar(type), nchar(type))) < 4) %>%
  select(-type)
annotation <- read_tsv('TAIR10_GFF3_genes.gff', col_names = c('chromosome', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes')) %>%
  filter(grepl('protein_coding_gene', attributes, fixed = TRUE)) %>%
  mutate(
    chromosome = sub('Chr', '', chromosome),
    gene = sub(';.*', '', substr(attributes, 4, 20))
  ) %>%
  select(gene, chromosome, strand)

convergent_gene <- sapply(seq_along(annotation$gene), function (x) ifelse(annotation$strand[x] == '+', annotation$gene[x + 1], annotation$gene[x - 1]))
names(convergent_gene) <- annotation$gene
### merge data with annotation ###
data_sense_annotated <- data_sense %>%
  inner_join(annotation, by = 'gene')

data_antisense_annotated <- data_antisense %>%
  inner_join(annotation, by = 'gene') %>%
  mutate(
    antisense_gene = gene,
    gene = convergent_gene[antisense_gene],
    antisense = TRUE
  ) %>%
  inner_join(annotation, by = 'gene', suffix = c('_antisense', '_sense')) %>%
  filter(chromosome_antisense == chromosome_sense & strand_antisense != strand_sense) %>%
  select('chromosome' = chromosome_sense, coord, gene, 'PAT' = `WT PAT#`, 'strand' = strand_sense, antisense)
### join sense and antisense data ###
data_both <- data_sense_annotated %>%
  mutate(
    antisense = FALSE
  )
colnames(data_both) <- c("coord","gene","PAT","chromosome","strand","antisense")
data_both <- data_both %>% bind_rows(data_antisense_annotated)



### keep only the most abundant PAC for each gene ###
# if there is a tie for the most abundant PAC, only the one yielding the shortest 3' UTR is kept
top_PACS <- data_both %>%
  group_by(chromosome, strand, gene) %>% slice_max(
    order_by = PAT,
    n = 1
  ) %>%
  mutate(
    helper = if_else(strand == '+', coord, -coord)
  ) %>%
  slice_min(
    order_by = helper,
    n = 1
  ) %>%
  ungroup() %>%
  group_by(across(-antisense)) %>%
  summarise(
    antisense = all(antisense)
  ) %>%
  select(-helper) %>% 
  ungroup()  %>% mutate(place = rep("1"))



second_top <- data_both %>%
  group_by(gene) %>% 
  arrange(by=desc(PAT),.by_group = TRUE) %>%
  slice_max(
    order_by = PAT,
    n = 2, with_ties=FALSE) %>% 
  mutate(percent_of = PAT/sum(PAT)*100) %>%
  mutate(max = max(PAT)) %>%
  mutate(helper = if_else(PAT == secondMax(PAT),TRUE,FALSE)) %>%
  filter(helper == TRUE) %>%
  filter(percent_of > 30) %>%
  select(-percent_of,-max,-helper) %>% 
  ungroup() %>%
  group_by(across(-antisense)) %>%
  summarise(
    antisense = all(antisense)
  ) %>%
  ungroup() %>% mutate(place = rep("2"))


combined <- rbind(top_PACS,second_top) %>% arrange(by=gene) %>% mutate(gene=if_else(place==2,paste(gene,"_2",sep = ''),gene)) 


combined %>%
  mutate(
    start = if_else(strand == '+', coord - 151, coord - 20),
    end = if_else(strand == '+', coord + 19, coord + 150),
    score = paste0(if_else(antisense, 'a', 's'), PAT)
  ) %>%
  select(chromosome, start, end, gene, score, strand) %>%
  arrange(chromosome, start, end) %>%
  write_tsv('Arabidopsis_PACs_Thomas_plus_second.bed', col_names = FALSE)


arab_thomas_2 <- read.table("Arabidopsis_PACs_Thomas_plus_second.bed", header=FALSE)
arab_thomas_2$V1<- sub("^", "Chr", arab_thomas_2$V1 )


Not_in_thomas <- arab_170 %>% filter(!(V4 %in% arab_thomas_2$V4))

Minus <- Not_in_thomas %>% filter(V6=="-") %>% mutate(V2=V2-1) %>% mutate(V3=V3-1) 
Plus <- Not_in_thomas %>% filter(V6=="+") 


mashed <- rbind(Minus,Plus)

Final_arab <- rbind(arab_thomas_2,mashed) %>% arrange(by=V1,V2,V4) 

duplicated_arab <- rbind(arab_thomas_2,mashed) %>% group_by(V1,V2,V3,V6) %>% filter(n()>1) %>% summarise(newName = paste(V4, collapse=";")) %>% mutate(score=rep(".")) %>% select(V1,V2,V3,newName,score,V6)
colnames(duplicated_arab) <- c("V1","V2","V3","V4","V5","V6")

to_remove <- rbind(arab_thomas_2,mashed) %>% group_by(V1,V2,V3,V6) %>% filter(n()>1)

tmp <- Final_arab %>% filter(!(V4 %in% to_remove$V4))

Final_arab_2 <- rbind(duplicated_arab,tmp) %>% arrange(by=V1,V2,V4) %>% write_tsv('Arabidopsis_PACs_Thomas_plus_TAIR.bed', col_names = FALSE)


