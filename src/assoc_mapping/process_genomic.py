# This script process the "genomic" output from populations,
# turning the genotype encoding into a proportion of homozygousity
# and removing SNPs that are homozygous in mothers
# from their respective offspring
# Cyril Matthey-Doret
# 11.08.2017

import numpy as np
import pandas as pd

# Missing bases are encoded as 0
# Homozygous genotypes are: 1,5,8,10
# Heterozygous genotypes are all others (except 0)

def gen_decode(row):
    """"
    This function decodes numeric genotypes and
    replaces it with E (heterozygous), O (homozygous)
    or M (missing).
    :param row: A pandas dataframe row, respecting the
    genomics output structure from STACKS populations module.
    :returns: A numpy array containing the genotype letters.
    """

    genodict = {}
    for code in range(11):
        # Building dictionary for numeric genotype translation
        if code in [1,5,8,10]:  # Homozygous codes
            genodict[code] = 'O'
        elif code == 0:  # Missing
            genodict[code] = 'M'
        else:
            genodict[code] = 'E'  # All others are heterozygous
    out_row = [genodict[i] for i in row[:]]  # Translating genotypes

    return np.array(out_row)


def mother_hom(geno, pop):
    """
    This function runs on a numpy array that has already been
    transformed with gen_decode and sets SNP that are homozygous in
    mothers to missing in their whole family.
    :param pop: a dataframe containing individual names and their respective
    families. The names need to be in the same order as the columns in geno.
    :param geno: a numpy array that will be processed
    """

    for f in np.unique(pop.Family):  # Iterate over mothers
        fam = pop.loc[pop.Family == f,:]  # Subssetting samples from family
        mother_idx = fam.index[fam.Generation=='F3'].tolist()  # Get mother idx
        fam_SNP = np.where(geno[mother_idx]=='O')[0]  # hom. mother SNPs
        # Change those SNPs to M in all indv with same family
        geno.loc[fam_SNP,fam.index] = 'M'

    return geno


def prop_hom(geno):
    """
    This function computes the proportion f homozygous individuals (columns) at
    each SNP (row) in a numpy array containing decoded allelic state (O,E,M).
    :param geno: Numpy array with SNPs as rows and individuals as columns.
    :returns: a new array with the proportion of homozygous individuals at each
    SNP and the number of individuals where it was present.
    """
    N = geno.shape[1]  # Number of samples
    # Counting how many individuals are used to compute proportion at each SNP
    n_prop = geno.apply(lambda r: N - sum(r.str.count('M')),axis=1)
    prop = geno.apply(lambda r: (float(sum(r.str.count('O')))/
                      max(sum(r.str.count(r'[OE]')),1)),axis=1)
    return pd.DataFrame({"N.Samples":n_prop,"Prop.Hom":prop})


### LOADING AND PROCESSING DATA ###

in_path = "../../data/populations/d-3_r-80/"  # argv1
indv_path = "../../data/individuals"  # argv2
out_path = "../../data/assoc_mapping/prop_hom_fixed_sites.tsv"  # argv3
genomic = pd.read_csv(in_path + "batch_0.genomic.tsv", sep='\t', header=None,
                      skiprows=1)
# Preparing data structure to match sample names and families with columns
indv = pd.read_csv(indv_path, sep='\t')  # Family and sex info
samples = pd.read_csv(in_path + "batch_0.sumstats.tsv",
            sep='\t',nrows=2, header=None)  # Names in correct order
# Concatenating populations
names = samples.iloc[:,1][0].split(',') + samples.iloc[:,1][1].split(',')
names = pd.DataFrame({'Name':names})
# Adding family and sex, keeping order
pop = names.merge(indv,on='Name',how='left')

### RUNNING CODE ###
genomic = genomic.iloc[:10,:]
gen_indv = genomic.iloc[:,3:].T.reset_index(drop=True).T  # only samples cols
# Decoding numeric genotypes into states (het, hom, missing)
state = gen_indv.apply(lambda r: gen_decode(r), axis=1)
clean = mother_hom(state, pop)  # Removing SNPs that are homozygous in mothers
prop = prop_hom(clean)

### SAVING OUTPUT ###
print prop.shape
print genomic.shape
#prop.loc[:,["Locus.ID","Chr","BP"]] = genomic.iloc[:,0:2]
prop = pd.concat([prop,pd.DataFrame(genomic.iloc[:,0:2],
                         columns=["Locus.ID","Chr","BP"],
                         index=prop.index)], axis=1)
prop.to_csv(out_path, sep='\t', index=False)
