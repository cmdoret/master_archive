# This script allows to differentiate between haploid and diploid males based on
# heterozygosity. It generates a list with each male and their ploidy.
# Cyril Matthey-Doret
# 11.05.2017

import pandas as pd
from sys import argv

# 1: read data/mother_fam.txt to extract families of mothers
vcf_sum = pd.read_csv(argv[1], sep='\t')
# Using command line argument to load VCF file summary
fam_sum = pd.read_csv("data/individuals",sep='\t')

fam_sum = fam_sum.merge(vcf_sum, left_on='Name', right_on='INDV', how='left')
fam_sum.dropna(how='any',axis=0,inplace=True)
females = fam_sum[fam_sum.Sex == 'F']

# 2: compute mean and standard deviation Fis of females in each family
Fis_ref = pd.DataFrame({'MEAN':females.groupby(['Family'])['F'].mean(),
                         'STD':females.groupby(['Family'])['F'].std()})

# 3: In each family: males with heterozygosity < mean(F) - N*stdev(F) -> haploid

males = fam_sum[fam_sum.Sex == 'M'].groupby(['Family'])

haplo = males.apply(lambda g:
                      g[g['F'] >= (Fis_ref.loc[g.name,'MEAN'] +
                                   2*Fis_ref.loc[g.name,'STD'])])


# 4: output list of males with ploidy information
out_haplo = fam_sum.loc[fam_sum.Sex == 'M',['Name','Family']]
out_haplo = out_haplo.merge(haplo,on='Name',how='inner')[['Name','Family_x']]
out_haplo.rename(columns={'Family_x':'Family'},inplace=True)
out_haplo.to_csv('data/haploid_males', sep='\t')
