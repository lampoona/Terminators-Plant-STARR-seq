#install.packages("dplyr")
library(dplyr)


arab <- read.table("TAIR10_GFF3_genes_three_prime_UTR.gff", header=FALSE)
arab$sort <- substr(arab$V4,5,9)
arab$sort <- as.numeric(arab$sort)
newdata <- arab[order(arab$V1,arab$sort),]
fake <- newdata %>% select(-sort)
#fake[duplicated(fake$V4) | duplicated(fake$V4, fromLast = TRUE),]
write.table(fake, file="TAIR10_GFF3_genes_three_prime_UTR.bed", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
