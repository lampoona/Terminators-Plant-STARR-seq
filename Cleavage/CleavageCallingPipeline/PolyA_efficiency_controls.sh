#!/bin/bash

#have to load a different module
module load samtools/1.10

Name_file=control_gene_names.txt
names=$(cat $Name_file)
for gene in $names
do
	#cd ${gene}
	samtools view -b -X processed_Light_plasmid_sorted_filtered.bam processed_Light_plasmid_sorted_filtered.bam.bai ${gene} > ${gene}/${gene}.bam
	samtools view ${gene}/${gene}.bam | cut -f1,3,10 | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' > ${gene}/${gene}_RN_UMI_GENE_SEQ.txt 
	variable=$(grep ${gene} gene_names_and_length.txt | awk '{print $2}') 
	awk -v gene_name=${gene} -v size=${variable} '{arr[$2]=length($4)} END {for (len in arr) {all++; if (arr[len] <= size) cut++} print gene_name,cut/all}' ${gene}/${gene}_RN_UMI_GENE_SEQ.txt >> polyEff_controls.txt
	rm ${gene}/${gene}.bam
done  
