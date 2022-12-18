import re,os,sys,glob
import mycustom

def reader(path):
        datafile=open(path,"r")
        data=datafile.readlines()
        datafile.close()
        Data=[]
        for i in data:
                Data+=[i.strip().split('\t')]
        return Data

AminoAcidsL=["G","A","L","M","F","W","K","Q","E","S","P","V","I","C","Y","H","R","N","D","T"]


name=sys.argv[1]
kmersS=set()
# folders = all proteomes to be studies
# name species to be excluded
folders=["Archaea_proteomes","Bacteria_proteomes","Viruses_proteomes","Eukaryota_proteomes"]
for one in folders:
	files=glob.glob(one+"/*_6mers")
	if name not in one:
		for two in files:
			DataL=reader(two)
			KmersL=[v[0] for v in DataL]
			for kmer in KmersL:
				kmersS.add(kmer)
# kmerS is the intersection of all unique kmer peptides from all species except the excluded species

filed=glob.glob("Bacteria_proteomes/"+name+"_*"+".fasta_6mers")[0]
# filed = one file containing the excluded species (name)

spp_kmersS=set()
for one in filed:
	for k in reader(one):
		spp_kmersS.add(k[0])

#spp_specific=list(spp_kmersS.difference.kmersS)
spp_specific=list(spp_kmersS.difference(kmersS))
datafile=open(name+"_Qprimes_6mers2.txt","w")
for k in spp_specific:
	datafile.write(k+'\n')
datafile.close()

# Calculate quasi-primes for given species
