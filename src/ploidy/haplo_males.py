# This script allows to differentiate between haploid and diploid males based on
# heterozygosity. It generates a list with each male and their ploidy.
# Cyril Matthey-Doret
# 11.05.2017

import pandas as pd
from sys import argv

vcf_sum = pd.read_csv(argv[1])
