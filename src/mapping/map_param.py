import matplotlib.pyplot as plt
import pandas as pd

"""
This script generates plots from custom generated summaries of the BWA
mapping output with different parameters. This is useful to visualize which
combinations of parameters yield the best results.
Cyril Matthey-Doret
02.04.2017
"""

# Importing mapping output summary table
map_sum = pd.read_csv("mapstats.csv",header=0)  # Loading mapping summaries

# Splitting by mapping algorithm
aln = map_sum.loc[map_sum['alg'] == 'aln']
mem = map_sum.loc[map_sum['alg'] == 'mem']
print(aln.loc["MM"])
pd.plot(aln["MM"],aln[])
