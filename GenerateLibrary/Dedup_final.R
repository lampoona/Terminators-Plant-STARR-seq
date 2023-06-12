
library(tidyverse)
library(tidyr)
library(dplyr)

all_sequences <- read.table("RR_Final_terminator_sequences.fa",header=FALSE) %>%  mutate(V1= gsub('>', '', V1))

duplicated_names  <- all_sequences %>% group_by(V1) %>% arrange(by=V1) %>% 
  filter(n()>1)

add_in <- duplicated_names %>% separate(V1, c("A", "B")) %>% mutate(place = if_else(row_number() %% 2 == 1,1,2)) %>% mutate(newname=paste(A,place,B,sep="_")) %>% select(newname,V2)
colnames(add_in) <- c("V1","V2")
add_in_1 <- add_in %>% mutate(V1= if_else(substr(V1,1,2)=="Zm",paste(substr(V1,1,16),"_CDS",sep=""),V1))
removed <- all_sequences %>% filter(!(V1 %in% duplicated_names$V1))

added <- rbind(removed,add_in_1)

added %>% write_delim('Final_terminator_sequences.fa', col_names=FALSE,delim="\n")
