library(tidyverse)
library(readr)
library(readxl)
library(dplyr)

Final_arab_2 <- read.table('Arabidopsis_PACs_Thomas_plus_TAIR.bed', header = FALSE)

Arab_sequences <- read.table("R_Arabidopsis_PACs_Thomas_plus_TAIR.fa",header=FALSE) %>%  mutate(V1= gsub('>', '', V1))

duplicated_seq <- Arab_sequences %>% group_by(V2) %>% arrange(by=V2) %>% 
  filter(n()>1) #%>%  summarise(distinct_gene_name = n_distinct(V1))

duplicated_f <- duplicated_seq %>% group_by(V2) %>% mutate(newName = paste(V1, collapse=";")) %>% select(V1,newName,V2) 
colnames(duplicated_f) = c('V4','V7','V8')

to_remove_again <- Final_arab_2 %>% filter(V4 %in% duplicated_f$V4)

to_add_back <- Final_arab_2 %>% filter(V4 %in% duplicated_f$V4) %>% left_join(duplicated_f,by='V4') %>% arrange(by=V8) %>% mutate(helper=row_number()) %>% filter(helper %% 2 == 1) %>% select(V1,V2,V3,V7,V5,V6) 
colnames(to_add_back) = c("V1","V2","V3","V4","V5","V6")

tmp <- Final_arab_2 %>% filter(!(V4 %in% to_remove_again$V4))

Final_arab_3 <- rbind(tmp,to_add_back) %>% arrange(by=V1,V2,V4) %>% write_tsv('Arabidopsis_PACs_Thomas_plus_TAIR.bed', col_names = FALSE)

