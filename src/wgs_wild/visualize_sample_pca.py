# This script performs PCA on the whole genome and on the different regions with
# high nucleotidic diversity to visualize which samples sample relatedness
# in the different regions.
# Cyril Matthey-Doret
# 09.04.2018

import numpy as np
import pandas as pd
from sklearn.decomposition import RandomizedPCA

snp_df = pd.read_table('data/wgs_wild/variant/ez.chr.wild.matrix.txt',
                       sep='\t', header=None)
