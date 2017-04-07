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
single_aln = aln.loc[:,"single"]/aln.loc[:,"tot"]
single_mem = mem.loc[:,"single"]/mem.loc[:,"tot"]

# plotting
"""
fig = plt.figure()
ax = fig.add_subplot(121)
ax.plot(aln.loc[:,"MM"],single_aln)
ax.set(title="aln single hits",
       ylabel="proportion of single mapped reads",
       xlabel="number of mismatches allowed")
bx = fig.add_subplot(122)
bx.plot(mem.loc[:,"K"],single_mem)
bx.set(title="mem single hits",
       ylabel="proportion of single mapped reads",
       xlabel="minimum seed length")
"""
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(aln.loc[:,"MM"],single_aln)
ax.set(title="aln single hits",
       ylabel="proportion of single mapped reads",
       xlabel="number of mismatches allowed")
plt.savefig("mapstats.pdf")
