#Sayeh Gorjifard January 2022
#Program takes a bed file and returns a bed file with only the last "exons" for each gene. It's critical to the program that the bedfile is resorted by chromosome and then gene name after it has already been sort-bed. 

import sys
import string
import os
import re
import pandas as pd
import csv 

file = sys.argv[1]

base=os.path.basename(file)
name=os.path.splitext(base)[0].split('_')[0]
print(base)
print(name)

df = pd.read_csv(file, sep='\t', comment='t', header=None)
header = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
df.columns = header[:len(df.columns)]

#initialize the first row
chrom = df.iloc[0,0]
name = df.iloc[0,3]
start = df.iloc[0,1]
end = df.iloc[0,2]
score=df.iloc[0,4]
strand = df.iloc[0,5]

final_list = []

df1 = df.iloc[1:df.shape[0],:]

#
for row in df1.itertuples():
	if row.strand == '+':
		if row.name == name:
			end = row.chromEnd
			start = row.chromStart
		else: 
			list1 = [chrom,start,end,name,score,strand]
			final_list.append(list1)
			chrom = row.chrom
			name = row.name
			start = row.chromStart
			end = row.chromEnd
			strand = row.strand
			score = row.score
	elif row.strand == '-':
		if row.name != name:
			list2 = [chrom,start,end,name,score,strand]
			final_list.append(list2)
			chrom = row.chrom
			name = row.name
			start = row.chromStart
			end = row.chromEnd
			strand = row.strand
			score = row.score
		
list2 = [chrom,start,end,name,score,strand]
final_list.append(list2)

data = pd.DataFrame(final_list, columns= header)
filename= str(name) +"_3UTR_last_exon.bed"
data.to_csv("TAIR10_3UTR_last_exon.bed", sep ='\t', quoting=csv.QUOTE_NONE,header=False, index=False)
