# Used to decide which haploid offspring to merge in merge_haplo.sh
# takes a pedigree file output by ParentCall2 from LEPMAP3 as input.
# outputs a list of pairs in the format
# p1k p2k p3k
# p1d p2d p3d
# e1 e2 e3
# Where pNk is the member of pairN that must be kept and pNd is the member that must be dropped
# Members are represented by 1-indexed column number
# Samples to be excluded (parents and unpaired offspring) are on the 3rd line (e1 e2...)
# CMD
# 18.06.2018

import pandas as pd
import numpy as np

# Read pedigree part of the file, ignoring auto-commented lines
ped = pd.read_csv("data/linkage_map/lib13/F4_geno.call", nrows=6, comment="#", header=None, sep='\t')
# Remove irelevant pedigree variables and coordinate infos
ped = ped.iloc[[0,1,3],2:]

# Transpose rows and cols
ped = ped.T
# Rename columns to meaningful names
ped.columns = ["Family", "Name", "Parent"]
# Assign a number representing sample rank in the family (same as in pedigree)
ped["fam_idx"] = ped.groupby("Family").cumcount()+1
# Drop all samples with odd ranks
ped["fate"] = ped.fam_idx.map(lambda x: "d" if x % 2 else "k")

# Exclude samples without parents (i.e. parents themselves)
ped.fate = np.where(ped['Parent'] == "0", 'e', ped.fate)

# Exclude samples if fam_idx is odd AND max among offspring
def exclude_odd(df):
	df.fate[(df.fam_idx == max(df.fam_idx[df.Parent != "0"]))  & (df.fam_idx % 2)] = 'e'
	return df

ped = ped.groupby("Family").apply(exclude_odd)


# Retrieve index, make 1-indexed and write to output file
with open("sample_pairs.txt", 'w') as outFile:
	for f in ['k','d','e']:
		for item in list(map(lambda x: int(x)+1, ped.index[ped['fate'] == f].tolist())):
			outFile.write("%s\t" % item)
		outFile.write("\n")
