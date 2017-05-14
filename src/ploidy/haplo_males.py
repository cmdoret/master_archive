# This script allows to differentiate between haploid and diploid males based on
# heterozygosity. It generates a list with each male and their ploidy.
# Cyril Matthey-Doret
# 11.05.2017

import pandas as pd
from sys import argv

# 1: read data/mother_fam.txt to extract families of mothers
vcf_sum = pd.read_csv(argv[1], sep='\t')
# Using command line argument to load VCF file summary
vcf_sum['Family'] = vcf_sum.INDV.str.slice(start=3,stop=4)
# Extracting family from individual names (will only work on F4 individuals)
motherfam = pd.read_csv("../../data/mother_fam.txt",sep='\t')

# vcf_sum.loc[vcf_sum.INDV.isin(motherfam.Parent),'Family'] = mother_fam['Family']
fam_sum = vcf_sum.merge(motherfam, left_on='INDV', right_on='Parent', how='outer')
vcf_sum.merge(motherfam, left_on='INDV', right_on='Parent', how='left')
pd.Series().isnull
# 2: compute mean Fis of females in each family
Fis_stat = fam_sum.groupby(['Family'])['F'].mean()
print(Fis_stat)
# 3: compute standard deviation of Fis in each family

# 4: In each family: males with heterozygosity < meanF - NstdevF -> haploid
# 5: output list of males with ploidy information
