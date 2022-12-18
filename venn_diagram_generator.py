import re,os,sys,glob,itertools,pdb
import mycustom 
from Bio.Seq import Seq 
from collections import Counter
import pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt 
from scipy.stats import pearsonr,binom_test

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data


# Generate Venn diagrams
# Mention in github: Taxonomic quasi-primes can be generated in the same way
groups=["archaea_all_7mers.txt","eukaryota_all_7mers_from_parts.txt","viruses_all_7mers.txt","bacteria_all_7mers_from_parts.txt"]
for one in groups:
	kmers1L=[k[0] for k in reader(one)]
	print(one,len(kmers1L))
	for two in groups:
		if one!=two:
			kmers2L=[k2[0] for k2 in reader(two)]
			print(two,len(kmers2L))
			common=len(set(kmers1L).intersection(set(kmers2L)))
			print(one,two,common)
		
			for three in groups:
				if one!=three and two!=three:
					kmers3L=[k2[0] for k2 in reader(three)]
					print(three,len(kmers3L))
					common=len(set(kmers3L).intersection(set(kmers1L).intersection(set(kmers2L))))
					print(one,two,three,common)
			
					for four in groups:
						if one!=four and two!=four and four!=three:
							kmers4L=[k2[0] for k2 in reader(four)]	
							print(four,len(kmers4L))
							common=len(set(kmers4L).intersection(set(kmers3L).intersection(set(kmers1L).intersection(set(kmers2L)))))
							print(one,two,three,four,common)
