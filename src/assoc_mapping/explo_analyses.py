
# This script is used for preliminary analyses before association mapping.
# It computes the number of alleles which fit the CSD pattern in all
# individuals.
# Cyril Matthey-Doret
# 10.06.2017

import pandas as pd
import numpy as np
from sys import argv
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

## Loading data
#in_geno = argv[1]
in_geno = "../../data/populations/d-25_r-75/batch_0.haplotypes.tsv"
geno = pd.read_csv(in_geno,sep='\t')  # importing genotypes into df
consensus = re.compile(r'consensus|-')  # faster to compile regex first

# Excluding all rows with only consensus or missing
consensus_rows  = geno.iloc[:,2:].apply(
       lambda row : all([ consensus.match(e) for e in row ]), axis=1)
geno = geno[~consensus_rows]

traits = pd.read_csv('../../data/ploidy/thresholds/m2',sep='\t')
# Importing traits data (sex, ploidy, family...)

## Splitting data

# Subsetting males and females separately (excluding individuals)
# Note: comp for comp(lementary)
comp_mal = traits.loc[(traits.Sex=='F') | (traits.Ploidy == 'H'),"Name"]
comp_fem = traits.loc[(traits.Sex=='M') | (traits.Ploidy == 'H'),"Name"]

males = geno.drop(comp_mal,axis=1)  # Diploid males only
females = geno.drop(comp_fem,axis=1)  # Females only

# Initiating summary dataframes
males_hom_sum = males.iloc[:,0:2]  # Summary of homozygosity in males
females_het_sum = females.iloc[:,0:2]  # Summary of heterozygosity in females

"""
# Only SNPs always homozygous (males) or only heterozygous (females)
females_het_sum["het_full"] = females.iloc[:,2:].apply(
    lambda r: r.str.match(r'^[ACTG]+/[ACTG]+$').all(),axis=1)
males_hom_sum["hom_full"] = males.iloc[:,2:].apply(
    lambda r: r.str.match(r'^[ACTG]+$').all(),axis=1)
"""
## Computing statistics

# Extracting proportion of het./hom in fem/males for each SNP
females_het_sum["prop_femhet"] = females.iloc[:,2:].apply(
    lambda r: r.str.count(r'^[ACTG]+/[ACTG]+$'),axis=1).sum(axis=1)
males_hom_sum["prop_malhom"] = males.iloc[:,2:].apply(
    lambda r: r.str.count(r'^[ACTG]+$'),axis=1).sum(axis=1)

## Merging data

# Putting males and females info back together
csd_snp = females_het_sum.merge(males_hom_sum, on="Catalog ID", how="inner")
# Merging by SNP ID
csd_snp.drop("Cnt_x", 1, inplace=True)
csd_snp.rename(columns={"Cnt_y": "Cnt"},  inplace=True)
csd_snp["prop_CSD"] = (csd_snp.prop_femhet+csd_snp.prop_malhom)/2
# csd_snp.sort_values(by="prop_CSD", ascending=False)

## Visualizing data

fig, axScatter = plt.subplots(figsize=(5.5, 5.5))

# the scatter plot:
axScatter.scatter(csd_snp.prop_femhet, csd_snp.prop_malhom)
axScatter.set_aspect(1.)

# create new axes on the right and on the top of the current axes
# The first argument of the new_vertical(new_horizontal) method is
# the height (width) of the axes to be created in inches.
divider = make_axes_locatable(axScatter)
axHistx = divider.append_axes("top", 1.2, pad=0.1, sharex=axScatter)
axHisty = divider.append_axes("right", 1.2, pad=0.1, sharey=axScatter)

# make some labels invisible
plt.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
         visible=False)

# now determine nice limits by hand:
binwidth = 0.1
xymax = np.max([np.max(np.fabs(csd_snp.prop_femhet)),
                np.max(np.fabs(csd_snp.prop_malhom))])
lim = (int(xymax/binwidth) + 1)*binwidth

bins = np.arange(0, lim, binwidth)
axHistx.hist(csd_snp.prop_femhet, bins=bins)
axHisty.hist(csd_snp.prop_malhom, bins=bins, orientation='horizontal')

# the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
# thus there is no need to manually adjust the xlim and ylim of these
# axis.

#axHistx.axis["bottom"].major_ticklabels.set_visible(False)
for tl in axHistx.get_xticklabels():
    tl.set_visible(False)
#axHistx.set_yticks([0, 50, 100])

#axHisty.axis["left"].major_ticklabels.set_visible(False)
for tl in axHisty.get_yticklabels():
    tl.set_visible(False)
#axHisty.set_xticks([0, 50, 100])

plt.draw()
plt.show()
