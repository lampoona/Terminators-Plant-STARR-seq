#!/bin/bash

#have to load a different module
module load samtools/1.10

Name_file=gene_names_greaterthan20.txt
names=$(cat $Name_file)
for gene in $names
do
	#cd ${gene}
	samtools view -b -X processed_Light_plasmid_sorted_filtered.bam processed_Light_plasmid_sorted_filtered.bam.bai ${gene} > ${gene}/${gene}.bam
	samtools view ${gene}/${gene}.bam | cut -f1,3,10 | awk 'BEGIN{FS=OFS="\t"}{gsub("_","\t",$1)}1' > ${gene}/${gene}_RN_UMI_GENE_SEQ.txt 
	if [ -f "${gene}/${gene}_RN_UMI_GENE_SEQ.txt" ]; then
        echo "Processing ${gene}/${gene}_RN_UMI_GENE_SEQ.txt"
        
        # Add debugging statement to print the values of gene and size
        echo "Gene: ${gene}, Size: ${variable}"
      	variable=$(grep -w ${gene} gene_names_and_length.txt | awk '{print $2}')  
        # Your existing awk command here
        awk -v gene_name=${gene} -v size=${variable} '{arr[$2]=length($4)} END {for (len in arr) {all++; if (arr[len] <= size) cut++} print gene_name,cut/all}' ${gene}/${gene}_RN_UMI_GENE_SEQ.txt >> polyEff_revised.txt
        echo "Processing done for ${gene}"
    else
        echo "Error: File ${gene}/${gene}_RN_UMI_GENE_SEQ.txt not found."
    fi
	rm ${gene}/${gene}.bam
	rm ${gene}/${gene}_RN_UMI_GENE_SEQ.txt
done  

