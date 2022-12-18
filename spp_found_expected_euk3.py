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
        return data

PEPTIDE_LENGTH = 7
AminoAcidsL=["G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T"]
all_kmers = [''.join(p) for p in itertools.product(AminoAcidsL, repeat=PEPTIDE_LENGTH)]
Dic={};Dic2={};
for kmer in all_kmers:
	Dic[kmer]=0
	Dic2[kmer]=0

files=glob.glob("Eukaryota_proteomes/*7mers")
for one in files:
	print(one)
	DataL=reader(one)
	for k in DataL:
		try:
			Dic[k.strip()]+=1
		except:
			pass

# both "files" should have the same length
files=glob.glob("Sims_Eukaryota_proteomes/*7mers")
for one in files:
	print(one)
        DataL=reader(one)
        for k in DataL:
                try:
                        Dic2[k.strip()]+=1
                except:
                        pass


#total=len(glob.glob("Sims_Eukaryota_proteomes/*7mers"))
num_of_species=len(files)
perc_real=[];perc_simulated=[];
datafile=open("tables_number_of_species/Eukaryota_number_of_species_7mers_exp_obs__NEW.txt","w")
size_of_kmer_space = len(AminoAcidsL) ^ PEPTIDE_LENGTH
for k in Dic.keys():
	perc_real.append(Dic[k]/float(num_of_species))
	perc_simulated.append(Dic2[k]/float(num_of_species))
	if float(Dic[k]+Dic2[k])!=0:
		datafile.write(k+'\t'+str(Dic[k])+'\t'+str(Dic2[k])+'\t'+str(Dic[k]/float(Dic[k]+Dic2[k]))+'\t'+str(min(1,binom_test(Dic[k],Dic[k]+Dic2[k],0.5)*size_of_kmer_space))+'\n')
	else:
		datafile.write(k+'\t'+str(Dic[k])+'\t'+str(Dic2[k])+'\t'+str(0)+'\t'+str(min(1,binom_test(Dic[k],Dic[k]+Dic2[k],0.5)*size_of_kmer_space))+'\n')
datafile.close()
