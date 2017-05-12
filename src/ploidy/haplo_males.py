# This script allows to differentiate between haploid and diploid males based on
# heterozygosity. It generates a list with each male and their ploidy.
# Cyril Matthey-Doret
# 11.05.2017

import pandas as pd
from sys import argv

vcf_sum = pd.read_csv(argv[1], sep='\t')
motherfam = pd.read_csv("../../data/mother_fam.txt",sep='\t')
fam_sum = pd.merge(vcf_sum, motherfam, how='outer', left_on="INDV", right_on="Parent",
         left_index=False, right_index=False, sort=False,
         copy=True, indicator=False)
print(fam_sum)
# 1: read data/mother_fam.txt to extract families of mothers
# 2: compute mean Fis of females in each family
# 3: compute standard deviation of Fis in each family
# 4: In each family: males with heterozygosity < meanF - NstdevF -> haploid
# 5: output list of males with ploidy information
