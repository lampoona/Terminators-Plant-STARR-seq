#Sayeh Gorjifard January 2022
#File takes a fasta file with a given range (must be known) and output an oligo count for each position in the range


import sys
import string
import re
import pandas as pd
import csv 
import numpy as np
import os

file = sys.argv[1]

base=os.path.basename(file)

size = int(sys.argv[2])

def createList(r1, r2):
	return list(range(r1, r2+1))

with open (file, 'r') as f:
	#give range of positions where 1 is the -500 and 1000 is the +500 position relative to 3UTR

	A = [0] * size
	T = [0] * size
	G = [0] * size
	C = [0] * size

	total_sequences= 0

	flines = f.readlines()
	for fline in flines:
		if fline[0:1] == ">":
			total_sequences += 1
			continue
		else: 
			nline = fline
			for x in range(size):
				if nline[x] == "A":
					A[x] += 1
				if nline[x] == "T":
					T[x] += 1
				if nline[x] == "G":
					G[x] += 1
				if nline[x] == "C":
					C[x] += 1

r1, r2 = 1, size
data = {'positions' : createList(r1,r2), 'A_count': A, 'T_count': T, 'G_count': G, 'C_count': C}

# Create DataFrame
df = pd.DataFrame(data)

dfname=os.path.splitext(base)[0]+"_oligo_frequency_table.txt"

df.to_csv(dfname, sep ='\t', quoting=csv.QUOTE_NONE, header=True, index=False)







