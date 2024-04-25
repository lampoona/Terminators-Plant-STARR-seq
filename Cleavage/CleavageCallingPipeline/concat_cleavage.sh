#!/bin/bash

Name_file=gene_names_greaterthan20.txt
names=$(cat $Name_file)

for gene in $names
do 
	awk -v r=$gene 'OFS="\t" {if (NR>1) print $1,$2,$3=r}' ${gene}/${gene}_PACs.tsv >> cleavage_all.tsv
done  

