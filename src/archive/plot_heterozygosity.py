
# This script uses the output sumstats file from the populations module of the
# STACKS suite to exclude haploid males based on computed heterozygosity levels.

# Cyril Matthey-Doret
# 30.04.2017

import pandas as pd
import re
import matplotlib.pyplot as plt

vcf = pd.read_csv('../../data/populations/d-25_r-75/batch_0.vcf',sep='\t',
                  skiprows=9)

vcf.iloc[:,9:] = vcf.iloc[:,9:].apply(lambda x: x.str.split(':').str.get(0))

hom_freq = vcf.iloc[:,9:].apply(lambda x: x.str.match(r'([^\.])/\1').sum())
het_freq = vcf.iloc[:,9:].apply(lambda x: x.str.match(r'([^\.])/[^\1\.]').sum())

heterozygosity = het_freq/hom_freq
heterozygosity.to_csv('heterozygosity',sep="\t")
"""
haplotable = pd.read_csv("../../data/sstacks/batch_0.haplotypes.tsv",sep="\t")

haplotable = haplotable.drop(["Catalog ID","Cnt"],axis="columns")
haplosum = haplotable.loc[:].applymap(lambda x: str.count(x,'/')).mean()
haplosum.sort_values(inplace=True)
haplomales = haplosum[haplosum.index.str.contains(pat='_')]
tmp_dict = {"samples":haplomales.index,"het":haplomales.values}
df1 = pd.DataFrame(tmp_dict)
df1 = df1.assign(family=haplomales.index.str.slice(start=3,stop=4))
df1.set_index(['samples'],inplace=True)
df1 = df1.groupby('family')


rowlength = int(df1.ngroups/2)                         # fix up if odd number of groups
fig, axs = plt.subplots(figsize=(9,4),
                        nrows=2, ncols=rowlength,     # fix as above
                        gridspec_kw=dict(hspace=0.4)) # Much control of gridspec

targets = zip(sorted(df1.groups.keys()), axs.flatten())

for i, (key, ax) in enumerate(targets):
    ax.bar(range(len(df1.get_group(key).values[:,0])),
           df1.get_group(key).values[:,0])
    ax.set_xticklabels(df1.get_group(key).index.str.slice(start=3,stop=8)
                       ,rotation='vertical')
    ax.set_title('Family %s'%key)
ax.legend()
plt.rcParams["figure.figsize"] = (20,20)
plt.savefig("heterozygosity.pdf")
plt.show()
"""
