# This script process the "genomic" output from populations,
# turning the genotype encoding into a proportion of homozygosity
# and removing SNPs that are homozygous in mothers
# from their respective offspring. Since genomics output file from
# populations are typically huge, the script takes quite a long time
# to run.
# Cyril Matthey-Doret
# 11.08.2017

import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count  # Parallel computing support
from functools import partial  # "freeze" arguments when mapping function
from sys import argv

# Genotype encoding in genomic output from populations:
# Missing bases are encoded as 0
# Homozygous genotypes are: 1,5,8,10
# Heterozygous genotypes are all others (except 0)

def gen_decode(encoded):
    """"
    This function decodes numeric genotypes and
    replaces it with E (heterozygous), O (homozygous)
    or M (missing).
    :param encoded: A pandas dataframe, respecting the genomics output structure
    from STACKS populations module, only the genotype columns must be included.
    :returns: A dataframe containing the genotype letters.
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
    genodict[1090670464] = 'M'
    # Rarely, rows are filled this value. I assume this is a STACKS issue.
    decoded = encoded.apply(lambda r: np.array([genodict[i] for i in r])
                             , axis=1)
    return decoded


def mother_hom(geno, pop):
    """
    This function runs on a numpy array that has already been
    transformed with gen_decode and sets SNP that are homozygous in
    mothers to missing in their whole family. If the mother is not available,
    SNPs that are homozygous or missing in all offspring in the family are used
    instead as a proxy.
    :param pop: a dataframe containing individual names and their respective
    families. The names need to be in the same order as the columns in geno.
    :param geno: a numpy array that will be processed
    """
    f = np.unique(pop.Family)[1]
    for f in np.unique(pop.Family):  # Iterate over mothers
        fam = pop.loc[pop.Family == f,:]  # Subssetting samples from family
        mother_idx = fam.index[fam.Generation=='F3'].tolist()  # Get mother idx
        fam_SNP = np.where(geno[mother_idx]=='O')[0]  # hom. mother SNPs
        if not mother_idx:  # If the mother is not available
        # Use SNPs where no individual in the family is heterozygous instead
            fam_SNP = np.where(np.all(geno[fam.index].isin(['O','M']),
                                            axis=1))[0]
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
    # Computing proportion of homozygous individuals at each SNP (O/O+E)
    prop = geno.apply(lambda r: (float(sum(r.str.count('O')))/
                      max(sum(r.str.count(r'[OE]')),1)),axis=1)
    return pd.DataFrame({"N.Samples":n_prop,"Prop.Hom":prop})

def parallel_func(f, df, f_args=[], chunk_size=100):
    """
    Parallelizes a function that runs on a dataframe by splitting the dataframe
    into small chunks by rows and distributing chunks across several processes.
    :param f: Target function that will be parallelized
    :param df: pandas dataframe to be used as input
    :param f_args: optional arguments for the function to be parallelized. Need
    to be an iterable (list or tuple).
    :param chunk_size: size of the chunks in which df is split. Default=100
    :returns: the processed dataframe reconstructed by combining output from all
    processes
    """
    # Create pool of processes, size depends on number of core available
    pool = Pool(processes = cpu_count())
    tot_rows = df.shape[0]
    chunks = range(0,tot_rows, chunk_size)  # Start positions of chunks
    # Split df into chunks
    chunked_df = [df.iloc[c:(c+min(chunk_size,tot_rows)),] for c in chunks]
    func = partial(f, *f_args)  # Unpacking optional fixed arguments.
    result = pool.map(func, chunked_df)  # Mapping function to chunks.
     # Concatenating into single df. Order is preserved
    return pd.concat(result)


#========== LOADING AND PROCESSING DATA ==========#

in_path = argv[1]  # STACKS populations folder
out_path = argv[2] + "/prop_hom_fixed_sites.tsv"  # Path of output file
indv_path = "data/individuals"  # family and sex information
genomic = pd.read_csv(in_path + "/batch_0.genomic.tsv", sep='\t', header=None,
                      skiprows=1)
# Preparing data structure to match sample names and families with columns
indv = pd.read_csv(indv_path, sep='\t')  # Family and sex info
samples = pd.read_csv(in_path + "/batch_0.sumstats.tsv",
            sep='\t',nrows=2, header=None)  # Names in correct order
# Concatenating populations
names = samples.iloc[:,1][0].split(',') + samples.iloc[:,1][1].split(',')
names = pd.DataFrame({'Name':names})
# Adding family and sex, keeping order
pop = names.merge(indv,on='Name',how='left')

#========== RUNNING CODE ==========#

# genomic = genomic.iloc[:100,:]
gen_indv = genomic.iloc[:,3:].T.reset_index(drop=True).T  # only samples cols
# Decoding numeric genotypes into states (het, hom, missing)
state = parallel_func(gen_decode, gen_indv)
clean = mother_hom(state, pop)  # Removing SNPs that are homozygous in mothers
# Computing proportion of homozygous indv at each site
prop = parallel_func(prop_hom, clean)

#========== SAVING OUTPUT ==========#

prop = genomic.iloc[:,0:3].merge(prop,left_index=True, right_index=True)
prop.rename(columns={0:"Locus.ID",1:"Chr",2:"BP"}, inplace=True)
prop.to_csv(out_path, sep='\t', index=False)
