
# This script is used for preliminary analyses before association mapping.
# It computes descriptive statistics and allows to visualize how SNPS fit the
# CSD pattern in all individuals with a particular threshold. It takes 2 CL
# arguments: the haplotypes.tsv file from stacks populations and the threshold
# file containing detailed infos on individuals.
# Cyril Matthey-Doret
# 10.06.2017

import pandas as pd
import numpy as np
from sys import argv
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

################
# Loading data #
################

in_geno = argv[1]
thresh = argv[2]
# thresh = "../../data/ploidy/thresholds/m2"
# in_geno = "../../data/populations/d-25_r-75/batch_0.haplotypes.tsv"
geno = pd.read_csv(in_geno,sep='\t')  # importing genotypes into df
consensus = re.compile(r'consensus|-')  # faster to compile regex first

# Excluding all rows with only consensus or missing
consensus_rows  = geno.iloc[:,2:].apply(
       lambda row : all([ consensus.match(e) for e in row ]), axis=1)
geno = geno[~consensus_rows]

traits = pd.read_csv(thresh,sep='\t')
# Importing traits data (sex, ploidy, family...)

##################
# Splitting data #
##################

# Subsetting males and females separately (excluding individuals)
# Note: comp for comp(lementary)
comp_mal = traits.loc[(traits.Sex=='F') | (traits.Ploidy == 'H'),"Name"]
comp_fem = traits.loc[(traits.Sex=='M') | (traits.Ploidy == 'H'),"Name"]

males = geno.drop(comp_mal,axis=1)  # Diploid males only
females = geno.drop(comp_fem,axis=1)  # Females only

# Initiating summary dataframes
males_het_sum = males.iloc[:,0:2]  # Summary of heterozygosity in males
females_het_sum = females.iloc[:,0:2]  # Summary of heterozygosity in females

"""
# Only SNPs always homozygous (males) or only heterozygous (females)
females_het_sum["het_full"] = females.iloc[:,2:].apply(
    lambda r: r.str.match(r'^[ACTG]+/[ACTG]+$').all(),axis=1)
males_hom_sum["hom_full"] = males.iloc[:,2:].apply(
    lambda r: r.str.match(r'^[ACTG]+$').all(),axis=1)
"""

########################
# Computing statistics #
########################

# Extracting proportion of het. in fem/males for each SNP

# Number of homozygous SNPs in females
femhom_geno = females.iloc[:,2:].apply(
    lambda r: r.str.count(r'^[ACTG]+$'),axis=1).sum(axis=1)
# Number of heterozygous SNPs in females
femhet_geno = females.iloc[:,2:].apply(
    lambda r: r.str.count(r'^[ACTG]+/[ACTG]+$'),axis=1).sum(axis=1)

# Proportion calculation (het/(het+hom)). Done by dividing het count by
# total number of columns minus consensus and - .
females_het_sum["prop_femhet"] = femhet_geno / (femhet_geno+femhom_geno)

# Same operation in males
malhet_geno = males.iloc[:,2:].apply(
    lambda r: r.str.count(r'^[ACTG]+/[ACTG]+$'),axis=1).sum(axis=1)
malhom_geno = males.iloc[:,2:].apply(
    lambda r: r.str.count(r'^[ACTG]+$'),axis=1).sum(axis=1)

males_het_sum["prop_malhet"] = malhet_geno / (malhet_geno+malhom_geno)

################
# Merging data #
################

# Putting males and females info back together
csd_snp = females_het_sum.merge(males_het_sum, on="Catalog ID", how="inner")
# Merging by SNP ID
csd_snp.drop("Cnt_x", 1, inplace=True)
csd_snp.rename(columns={"Cnt_y": "Cnt"},  inplace=True)
csd_snp["prop_CSD"] = (csd_snp.prop_femhet+(1-csd_snp.prop_malhet))/2
# csd_snp.sort_values(by="prop_CSD", ascending=False)

####################
# Visualizing data #
####################
fe = csd_snp.prop_femhet; me = csd_snp.prop_malhet

fig, axScatter = plt.subplots(figsize=(5.5, 5.5))

# the scatter plot:
axScatter.scatter(fe, me,
                  c=fe + (1-me),
                  s=15, cmap='plasma')
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
binwidth = 0.05
xymax = np.max([np.max(np.fabs(fe)),
                np.max(np.fabs(me))])
lim = (int(xymax/binwidth) + 1)*binwidth

bins = np.arange(0, lim, binwidth)
axHistx.hist(fe, bins=bins)
axHisty.hist(me, bins=bins, orientation='horizontal')

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

# Adding labels to axes
axHisty.set_xlabel("Male SNPs")
axHistx.set_ylabel("Female SNPs")
axScatter.set_xlabel("Female heterozygosity")
axScatter.set_ylabel("Male heterozygosity")

# Displaying data summary on the plot
m1n = len(traits.Name[traits.Ploidy=='H'])
textstr = 'Threshold={0}\nM1N={1}\nM2N={2}\nF={3}\nSNPs={4}'.format(
    thresh.split('/')[-1], m1n, males.shape[1]-2, females.shape[1]-2,
    csd_snp.shape[0])
plt.text(1.3, 1.2, textstr, fontsize=9)
plt.draw()
plt.savefig("reports/lab_book/assoc_explo/" + thresh.split('/')[-1] + ".pdf")
# plt.show()
