#!/bin/bash

#script calls major PAC sites from the processed bam file. 
#have to load a different module
module load samtools/1.10

Name_file=gene_names_greaterthan20.txt
names=$(cat $Name_file)
for gene in $names
do
	mkdir ${gene}
	samtools view -b -X processed_Light_plasmid_sorted_filtered.bam processed_Light_plasmid_sorted_filtered.bam.bai ${gene} > ${gene}/${gene}.bam
	samtools view ${gene}/${gene}.bam | cut -f1,3,10 | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' > ${gene}/${gene}_RN_UMI_GENE_SEQ.txt 
	readnum=$(wc -l ${gene}/${gene}_RN_UMI_GENE_SEQ.txt | awk '{print $1}')
	if [ $readnum -ge 1000 ]; then Rscript PAC_Caller_2.R ${gene}/${gene}_RN_UMI_GENE_SEQ.txt; else Rscript revised_PAC_Caller_2.R ${gene}/${gene}_RN_UMI_GENE_SEQ.txt; fi
	rm ${gene}/${gene}.bam
	rm ${gene}/${gene}_RN_UMI_GENE_SEQ.txt
	rm Rplots.pdf
	rm *.sh.e*
	rm *.sh.o* 
done 
