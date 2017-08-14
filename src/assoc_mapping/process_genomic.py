# This script process the "genomic" output from populations,
# turning the genotype encoding into a proportion of homozygosity
# and removing SNPs that are homozygous or missing in mothers
# from their respective offspring. Since genomics output file from
# populations are typically huge, the script takes quite a long time
# to run although it is parallelized and will automatically use all
# core available on the host. It takes 2 arguments: the path to the
# folder containing populations files, used as input, and the folder where
# the output will be written.
# TODO: add sex-specific proportion of homozygosity and number of indv
# in prop_hom()
# TODO: make script compatible with per-family populations run.
# Cyril Matthey-Doret
# 11.08.2017

from os import path, walk
import argparse
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count  # Parallel computing support
from functools import partial  # "freeze" arguments when mapping function

#========== PARSING COMMAND LINE ARGUMENTS ==========#

# Genotype encoding in genomic output from populations:
# Missing bases are encoded as 0
# Homozygous genotypes are: 1,5,8,10
# Heterozygous genotypes are all others (except 0)
parser = argparse.ArgumentParser(description="This script processes the \
                                 'genomic' output from STACKS's populations \
                                 module to compute the proportion of homozygous\
                                 individuals at each genomic position.")
parser.add_argument('pop_files', type=str,
                                 help='Path to the folder containing \
                                 populations output files. Used as input',)
parser.add_argument('out', type=str,
                                 help='Folder where output will be written',)
parser.add_argument('--keep_all', action='store_true',
                                 help='Keep all SNPs, even if missing or\
                                 homozygous in the mother.',)

args = parser.parse_args()

#========== DEFINING FUNCTIONS ==========#

def unify_genomic(pop_path):
    """
    Reads populations "genomic" output files from each family's folder and
    gathers them into a single pandas DataFrame object.
    :param path: Path containing the family subfolders containing populations
    output files.
    :returns: A DataFrame containing all sites in all individuals.
    """

    # Adding trailing slash if not provided
    if pop_path[-1] != '/': pop_path += '/'
    first_fam = True
    #Iterating over family subfolders
    for subdir, dirs, files in walk(pop_path):
        if path.basename(subdir):
            # Read file from subfolder
            tmp = pd.read_csv(path.join(subdir, "batch_0.genomic.tsv"),
                              sep='\t', header=None, skiprows=1)
            if first_fam: # Assign dataframe on iteration run only
                united = tmp
            else:
                # Graft each family's samples onto final df as new columns
                # Using outer merge on chromosome, bp and Locus ID
                united = united.merge(tmp, how='outer', left_index=True,
                                      right_index=True, on=range(3))
            first_fam = False
    united = united.fillna(0)  # Change all NAs to the "missing" code
    return united

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
    decoded = encoded.apply(lambda r: np.array([genodict[i] for i in r]),axis=1)
    return decoded


def mother_hom(geno, pop):
    """
    This function runs on a numpy array that has already been
    transformed with gen_decode and sets sites that are homozygous/missing in
    mothers to missing in their whole family. If the mother is not available,
    sites that are homozygous or missing in all offspring in the family are used
    instead as a proxy.
    :param pop: a dataframe containing individual names and their respective
    families. The names need to be in the same order as the columns in geno.
    :param geno: a numpy array that will be processed
    """

    for f in np.unique(pop.Family):  # Iterate over mothers
        fam = pop.loc[pop.Family == f,:]  # Subssetting samples from family
        mother_idx = fam.index[fam.Generation=='F3'].tolist()  # Get mother idx
        fam_SNP = np.where(geno[mother_idx]!='E')[0]  # hom/missing mother sites
        if not mother_idx:  # If the mother is not available
        # Use sites where no individual in the family is heterozygous instead
            fam_SNP = np.where(np.all(geno[fam.index].isin(['O','M']),
                                            axis=1))[0]
        # Change those sites to M in all indv with same family
        geno.loc[fam_SNP,fam.index] = 'M'

    return geno


def prop_hom(pop, geno):
    """
    This function computes the proportion of homozygous individuals (columns) at
    each SNP (row) in a numpy array containing decoded allelic state (O,E,M).
    It computes this proportion both by sex, and in all individuals.
    :param geno: Numpy array with sites as rows and individuals as columns.
    :param pop: Dataframe containing the sex of each individual and its name.
    :returns: a Pandas DataFrame object with the proportion of homozygous
    females, males and all individuals at each site and the number of individuals
    where it was present.
    """

    # Number of males and females
    N = {sex:pop.Sex[pop.Sex == sex].shape[0] for sex in ['M','F']}
    # Get sample indices by sex
    idx = {sex:pop.index[pop.Sex == sex] for sex in ['M','F']}

    # Counting how many individuals are used to compute proportion at each SNP
    sample_size = {}  # Number of individuals in which each site is found
    hom = {}  # proportion of individuals in which each site is homozygous
    for sex in N:
        # Looping over sexes
        sample_size[sex] = geno.iloc[:,idx[sex]].apply(lambda r:
                                       N[sex] - sum(r.str.count('M')),axis=1)
        # Computing proportion of homozygous individuals at each SNP (O/O+E)
        hom[sex] = geno.iloc[:,idx[sex]].apply(lambda r:
                     (float(sum(r.str.count('O')))/
                      max(sum(r.str.count(r'[OE]')),1)),axis=1)

    # Building output dataframe with all relevant stats
    out_df = pd.DataFrame({
        "N.Samples": sample_size['F'] + sample_size['M'],
        "Prop.Hom": (sample_size['M'] * hom['M'] +
                    sample_size['F'] * hom['F']) /
                    (sample_size['F'] + sample_size['M']),
        "N.Males": sample_size['M'],
        "N.Females": sample_size['F'],
        "Prop.Hom.F": hom['F'],
        "Prop.Hom.M": hom['M']
        })
    return out_df

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
# Path to STACKS populations folder and output file
in_path = args.pop_files
out_path = args.out + "/prop_hom_fixed_sites.tsv"
indv_path = "data/individuals"  # family and sex information
genomic = pd.read_csv(path.join(in_path, "batch_0.genomic.tsv"),
                      sep='\t', header=None, skiprows=1)
# Preparing data structure to match sample names and families with columns
indv = pd.read_csv(indv_path, sep='\t')  # Family and sex info
samples = pd.read_csv(path.join(in_path, "batch_0.sumstats.tsv"),
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

# Will run unless user explicitly set the --keep_all parameter
if not args.keep_all:
    # Remove SNPs that are hom./missing in mothers from their family
    state = mother_hom(state, pop)
# Computing proportion of homozygous indv at each site
prop = parallel_func(prop_hom, state, f_args = (pop,))

#========== SAVING OUTPUT ==========#

# Merging Chromosomal positions with proportion of homozygosity into 1 df
prop = genomic.iloc[:,0:3].merge(prop,left_index=True, right_index=True)
prop.rename(columns={0:"Locus.ID",1:"Chr",2:"BP"}, inplace=True)
prop.to_csv(out_path, sep='\t', index=False)
