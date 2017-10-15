# This script uses EOM matrix made from the genomic output of populations to
# identify loci that are heterozygous in haploid males and stores them in a
# blacklist file.
# Cyril Matthey-Doret
# 13.10.2017

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Uses EOM matrix made from the \
                                 genomic output of populations to identify \
                                 loci that are heterozygous in haploid males \
                                 and stores them in a blacklist file.")

parser.add_argument('EOM_in', type=str, help='Path to the EOM genotype matrix.')
parser.add_argument('black_out', type=str,
                    help='Path to the output blacklist file.')
parser.add_argument('ploid', type=str, help='Path to the table with ploidy \
                    information.')
parser.add_argument('--haplo_het', type=int, help='Proportion (float between 0 \
                    and 1) of haploid males required to be heterozygous at a \
                    locus in order to exclude it.')
args = parser.parse_args()

# Load genotype matrix and ploidy table
EOM = pd.read_csv(args.EOM_in, sep='\t')
ploid = pd.read_csv(args.ploid, sep='\t')

# Subset haploid individuals
haplonames = ploid.loc[ploid.Ploidy=='H',"Name"]
haplo_EOM = EOM.loc[:,haplonames]


# Suppresses warning when numpy divides by 0
np.seterr(divide='ignore', invalid='ignore')

dff = {}
for t in ['O','M','E']:
    # Store arrays with number of E/O/M samples at each locus
    dff[t] = (haplo_EOM==t).T.sum().astype(float)

# Computing proportion of heterozygosity at each locus
het = np.divide(dff['E'], (dff['O'] + dff['E']))
het = het.fillna(0)

# Extract loci above heterozygosity threshold
het_idx = het.index[(het>haplo_het)]
blacklist = EOM.iloc[het_idx,0].unique()

# Saving blacklisted loci to file
pd.DataFrame(blacklist).to_csv(black_out, sep='\t', header=False, index=False)
