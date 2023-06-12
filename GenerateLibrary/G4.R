
library(tidyverse)
library(readr)
library(readxl)
library(dplyr)

#reading in G4s
Maize_G4 <- read.table("Maize_Top_PACs_G4.bed")
Arab_G4 <- read.table("Arabidopsis_PACs_Thomas_plus_TAIR_G4.bed")

#read in final sequence with RE removed and adapters added 
Arab_sequences_final <- read.table("RR_Arabidopsis_PACs_Thomas_plus_TAIR_final.fa",header=FALSE) %>%  mutate(V1= gsub('>', '', V1))
Arab_sequences_final$no_adapter <- str_sub(Arab_sequences_final$V2,16,185)

Maize_sequences_final <- read.table("RR_Maize_Top_PACs_Final_Final.fa",header=FALSE) %>% mutate(V1= gsub('>', '', V1))
Maize_sequences_final$no_adapter <- str_sub(Maize_sequences_final$V2,16,185)

#24 nt sequences of G4s 
AAAA <- "GGGATATGGGATATGGGATATGGG"
BBBB <- "CCCATATCCCATATCCCATATCCC"
ABBB <- "GGGATATGGGATATGGGATATCCC"
AABB <- "GGGATATGGGATATCCCATATCCC"

set.seed(1000)
Arab_G4_200 <- Arab_sequences_final %>% filter(!(V1 %in% Arab_G4$V1)) %>% slice_sample(n = 200)

tmp1 <- Arab_G4_200 %>%  mutate(take_out = str_sub(no_adapter,38,61))
tmp1$AAAA <- str_replace(tmp1$no_adapter, as.character(tmp1$take_out), AAAA)
tmp1$BBBB <- str_replace(tmp1$no_adapter, as.character(tmp1$take_out), BBBB)
tmp1$ABBB <- str_replace(tmp1$no_adapter, as.character(tmp1$take_out), ABBB)
tmp1$AABB <- str_replace(tmp1$no_adapter, as.character(tmp1$take_out), AABB)
tmp1$V1<-gsub("^",">",tmp1$V1)


RE.sites <- c('(G|C|^)GTCTC', '(G|^)AGAC(C|G)')
temp1 <- tmp1 %>% mutate(hasNo_RE_AAAA = !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp1$AAAA)) %>% mutate(hasno_RE_BBBB= !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp1$BBBB)) %>% mutate(hasno_RE_AABB = !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp1$AABB)) %>% mutate(hasno_RE_ABBB = !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp1$ABBB)) %>%
  mutate(Total = select(., hasNo_RE_AAAA:hasno_RE_ABBB) %>% rowSums(na.rm = TRUE)) %>% filter(Total == 4) %>%
  slice(1:100)


temp1 %>% select(V1,AAAA) %>% mutate(name = paste(V1,"AAAA",sep = "-")) %>% select(name,AAAA) %>% write_delim('Arab_AAAA_G4.fa', col_names = FALSE, delim = "\n")
temp1 %>% select(V1,BBBB) %>% mutate(name = paste(V1,"BBBB",sep = "-")) %>% select(name,BBBB) %>% write_delim('Arab_BBBB_G4.fa', col_names = FALSE, delim = "\n")
temp1 %>% select(V1,ABBB) %>% mutate(name = paste(V1,"ABBB",sep = "-")) %>% select(name,ABBB) %>% write_delim('Arab_ABBB_G4.fa', col_names = FALSE, delim = "\n")
temp1 %>% select(V1,AABB) %>% mutate(name = paste(V1,"AABB",sep = "-")) %>% select(name,AABB) %>% write_delim('Arab_AABB_G4.fa', col_names = FALSE, delim = "\n")


set.seed(1000)
Maize_G4_200 <- Maize_sequences_final %>% filter(!(V1 %in% Maize_G4$V1)) %>% slice_sample(n = 200)

#take out the 24 nucleotides around the 50 location
tmp2 <- Maize_G4_200 %>%  mutate(take_out = str_sub(no_adapter,38,61))
tmp2$AAAA <- str_replace(tmp2$no_adapter, as.character(tmp2$take_out), AAAA)
tmp2$BBBB <- str_replace(tmp2$no_adapter, as.character(tmp2$take_out), BBBB)
tmp2$ABBB <- str_replace(tmp2$no_adapter, as.character(tmp2$take_out), ABBB)
tmp2$AABB <- str_replace(tmp2$no_adapter, as.character(tmp2$take_out), AABB)
tmp2$V1<-gsub("^",">",tmp2$V1)

RE.sites <- c('(G|C|^)GTCTC', '(G|^)AGAC(C|G)')
temp2 <- tmp2 %>% mutate(hasNo_RE_AAAA = !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp2$AAAA)) %>% mutate(hasno_RE_BBBB= !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp2$BBBB)) %>% mutate(hasno_RE_AABB = !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp2$AABB)) %>% mutate(hasno_RE_ABBB = !grepl(paste0('(', paste0(RE.sites, collapse = ')|('), ')'), tmp2$ABBB)) %>%
  mutate(Total = select(., hasNo_RE_AAAA:hasno_RE_ABBB) %>% rowSums(na.rm = TRUE)) %>% filter(Total == 4) %>%
  slice(1:100)

temp2 %>% select(V1,AAAA) %>% mutate(name = paste(V1,"AAAA",sep = "-")) %>% select(name,AAAA) %>% write_delim('Maize_AAAA_G4.fa', col_names = FALSE, delim = "\n")
temp2 %>% select(V1,BBBB) %>% mutate(name = paste(V1,"BBBB",sep = "-")) %>% select(name,BBBB) %>% write_delim('Maize_BBBB_G4.fa', col_names = FALSE, delim = "\n")
temp2 %>% select(V1,ABBB) %>% mutate(name = paste(V1,"ABBB",sep = "-")) %>% select(name,ABBB) %>% write_delim('Maize_ABBB_G4.fa', col_names = FALSE, delim = "\n")
temp2 %>% select(V1,AABB) %>% mutate(name = paste(V1,"AABB",sep = "-")) %>% select(name,AABB) %>% write_delim('Maize_AABB_G4.fa', col_names = FALSE, delim = "\n")








