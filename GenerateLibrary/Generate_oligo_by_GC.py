#python script takes GC content, length of desired sequence, and how many sequences you want as parameters. out to whatever file 
import sys
import string
import re
import pandas as pd
import csv 
import numpy as np
import os
import random

random.seed(42)

GC_freq = float(sys.argv[1])

size = int(sys.argv[2])

how_many = int(sys.argv[3])

#generate GC content
#then generate the remaining nucleotides AT
#join them together and randomly shuffle
def random_seq_gc_per(length, gc_perc):
    nb_gc = round(gc_perc*length)
    gc_seq = random.choices("GC", k=nb_gc)
    at_seq = random.choices("AT", k=length-nb_gc)
    dna_seq = gc_seq + at_seq
    random.shuffle(dna_seq)
    res_dna_seq = ''.join(dna_seq)
    return res_dna_seq

# filter out sequences with a BsaI or Esp3I site (the cloning adapter adds a G in front of the sequence)
#RE.sites <- c('(G|C|^)GTCTC', '(G|^)AGAC(C|G)')

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


for i in range(how_many):
    seq = random_seq_gc_per(size, GC_freq)
    while(check_for_RE(seq) == True):
    	seq = random_seq_gc_per(size, GC_freq)
    print(">"+"GC_"+str(int(GC_freq*100))+"_"+ str(i+1))
    print(seq)


