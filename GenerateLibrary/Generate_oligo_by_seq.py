#Sayeh Gorjifard
#February 28 2022
#script to create random oligos that match the overal ACGT content of the FASTA fed to it 
#python script takes fasta sequence, length of desired sequence, and how many sequences you want as parameters, and name of output sequence array. Must out ">" to a file 
import sys
import string
import re
import pandas as pd
import csv 
import numpy as np
import os
import random

random.seed(42)

#first takes in FASTA
file = sys.argv[1]

#second takes in length of oligos
size = int(sys.argv[2])

#how many generated sequences you want
how_many = int(sys.argv[3])

#the name subheading in the output 
name = str(sys.argv[4])

def getsequence(file): 
	outlines = []
	with open (file, 'r') as f: 
		count = 0
		for line in f:
			count+=1
			if count %2 ==0:
				outlines.append(line.strip())
	sequence  = ''.join(outlines)
	return sequence 



def getNucleotide(sequence):
	length = len(sequence)
	
	A = sequence.count('A')
	C = sequence.count('C')
	G = sequence.count('G')
	T = sequence.count('T')

	A_freq = A/length
	C_freq = C/length
	G_freq = G/length
	T_freq = T/length

	ACGT = [A_freq,C_freq,G_freq,T_freq]
	return ACGT

#obtain ACGT nucleotide frequency list to create a shuffled version with same frequency
#then generate the remaining nucleotides AT
#join them together and randomly shuffle
def random_seq_gc_per(length, comp_list):
	A = comp_list[0]
	C = comp_list[1]
	G = comp_list[2]
	T = comp_list[3]
	#print(comp_list)
	nb_A = round(A*length)
	nb_C = round(C*length)
	nb_G = round(G*length)
	nb_T = 170-nb_A-nb_C-nb_G
	#print(nb_A+nb_C+nb_G+nb_T)
	A_seq = "A" * nb_A
	C_seq = "C" * nb_C
	G_seq = "G" * nb_G
	T_seq = "T" * nb_T
	combined = A_seq + C_seq + G_seq + T_seq
	l = list(combined)
	random.shuffle(l)
	result = ''.join(l)
	return result


# filter out sequences with a BsaI or Esp3I site (the cloning adapter adds a G in front of the sequence)
#RE.sites <- c('(G|C|^)GTCTC', '(G|^)AGAC(C|G)')
regex = r"(G|^)AGAC(C|G)"
regex2= r"(G|C|^)GTCTC"


def check_for_RE(sequence): 
	testBool = False 
	regex = r"(G|^)AGAC(C|G)"
	regex2= r"(G|C|^)GTCTC"
	Test1 = bool(re.search(regex, sequence))
	Test2 = bool(re.search(regex2, sequence))
	if (Test1 == True):
		testBool = True
	elif (Test2 == True):
		testBool = True
	else: 
		testBool = False
	return testBool


DNA_sequence = getsequence(file)
ACGT = getNucleotide(DNA_sequence)


for i in range(how_many):
    seq = random_seq_gc_per(size, ACGT)
    while (check_for_RE(seq) == True):
    	seq = random_seq_gc_per(size, ACGT)
    print(">"+name+"_ACGT_random_"+ str(i+1))
    print(seq)














