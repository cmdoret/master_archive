
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
from os import path, walk
import re
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

################
# Loading data #
################

in_geno = argv[1]
thresh = argv[2]
grouped = argv[3]
# thresh = "../../data/ploidy/thresholds/fixed"
# in_geno = "../../data/ploidy/vcftools/"

SNP_splitter = lambda x: pd.Series([i for i in x.split('_')])
# splits SNP IDs that have been compacted in <contig_bp> form

traits = pd.read_csv(thresh,sep='\t')
# Importing traits data (sex, ploidy, family...)
mother_hom = pd.DataFrame() # Used to store all SNPs homozygous in mothers
fam_geno = {}
mat_pat = re.compile(r'012')
# Filename pattern for genotype matrices produced by vcftools
for subdir, dirs, files in walk(in_geno):
    # Walking through file tree
    if path.basename(subdir):
        # Exclude folder itself, only taking subfolders
        mat_files = list(filter(mat_pat.search,files))
        pos_index = [i for i, s in enumerate(mat_files) if 'pos' in s][0]
        pos = pd.read_csv(path.join(subdir, mat_files.pop(pos_index)),
                          header=None, sep='\t', dtype='str')
        # reading file with positions of SNPs
        pos = pos.iloc[:,0].str.cat(pos.iloc[:,1], sep='_')
        # Concatenating contig ID and position on each row
        indv_index = [i for i, s in enumerate(mat_files) if 'indv' in s][0]
        indv = pd.read_csv(path.join(subdir, mat_files.pop(indv_index)),
                              header=None, sep='\t', squeeze=True)
        # reading file with individuals in family
        geno_mat = pd.read_csv(path.join(subdir, mat_files[0]), sep='\t',
                               header=None, names=pos,
                               usecols=range(1,len(pos)+1))
        geno_mat.index = indv.values
        # reading genotypes matrix. ignoring first column as it contains row
        # numbers. Also giving individuals names as row names and SNPs positions
        # as column names
        try:
            if grouped == 'F':
                mother = list(traits.Name[(traits.Generation=='F3') &
                             (traits.Family==path.basename(subdir))])
            else:
                mother = list(traits.Name[(traits.Generation=='F3')])
            for m in mother:
                mo_bool = geno_mat.loc[[m],:]%2==0
                # Boolean mask for SNPs that are homozygous in mother
                mother_fam = pd.DataFrame(mo_bool.columns[mo_bool.iloc[0,:]],
                columns=['merged'])
                mother_fam[['contig', 'bp']] = mother_fam['merged'].apply(
                    SNP_splitter)
                mother_fam['family'] = traits.loc[traits['Name']==
                                                  m,'Family'].iloc[0]
                mother_fam.drop('merged', axis=1, inplace=True)
                # Split IDs to extract contig and base pair infos
                mother_hom = mother_hom.append(mother_fam)
                # Append family blacklist to overall dataframe
                fam_samples = traits.loc[traits['Family']==mother_fam.family[0],
                                         'Name']

                geno_mat.loc[fam_samples,mo_bool.iloc[0,:]] = -1
                # excluding SNPs that are homozygous in mothers
        except IndexError:
            print("No SNP data for mother of family " + path.basename(subdir))
        fam_geno[path.basename(subdir)] = geno_mat
        # Storing genotype matrix of each family in the dictionary

mother_hom.sort_values(['family', 'contig', 'bp'], inplace=True)
mother_hom.to_csv('data/SNP_lists/' + path.basename(thresh) + '_hom_mother.txt',
                        index=False)

# Sorting SNPs by family and writing them to a text file

##################
# Splitting data #
##################

# Subsetting males and females separately (excluding individuals)
# Note: comp for comp(lementary)
males = {}
females = {}
for fam in sorted(fam_geno):  # Iterating over families
    if grouped == "F":
        fam_traits = traits.loc[traits.Family==fam]
    else:
        fam_traits = traits
    # Subsetting traits df for each family
    comp_mal = fam_traits.loc[(fam_traits.Sex=='F') |
                          (fam_traits.Ploidy == 'H'),"Name"]
    comp_fem = fam_traits.loc[(fam_traits.Sex=='M') |
                              (fam_traits.Ploidy == 'H'),"Name"]
    # Individuals to be excluded from diploid males and females respectively
    males[fam] = fam_geno[fam].drop(comp_mal,axis=0)
    # Genotypes of diploid males only
    females[fam] = fam_geno[fam].drop(comp_fem,axis=0)
    # Genotypes of females only

########################
# Computing statistics #
########################
SNP_sum = {}
for fam in fam_geno:
    SNP_sum[fam] = pd.DataFrame(columns=['Male het.', 'Female het.'],
                                index=fam_geno[fam].columns.values)
    # Initiating dataframe to store SNPs statistics
    SNP_sum[fam].loc[:,'Male het.'] = males[fam].apply(
        lambda c:float(len(c[(c > 0) & (c % 2)])) /
        float(max(list([1, len(c[c > 0])]))),axis=0)

    SNP_sum[fam].loc[:,'Female het.'] = females[fam].apply(
        lambda c:float(len(c[(c > 0) & (c % 2)])) /
        float(max(list([1,len(c[c > 0])]))),axis=0)
    # Summarizing each SNPs by proportion of males and females in which it
    # is heterozygous. Note: not including individuals where SNP is absent.
    # max(1,prop_het) used to prevent divisiion by zero.

###############
# Output data #
###############

CSD_full = pd.DataFrame()

for fam in fam_geno:
    SNP_sum[fam]["prop_CSD"] = (SNP_sum[fam]["Female het."] +
                                  (1-SNP_sum[fam]["Male het."]))
    # Prop. of males in which SNP is hom. + prop of females in which it is het.
    # Basically a measure of how CSD-like it is
    CSD_like = SNP_sum[fam][SNP_sum[fam]["prop_CSD"]==1].index
    # subset those segragating exactly like CSD
    CSD_fam = pd.DataFrame(CSD_like, columns=['merged'])
    # Store them in dataframe
    CSD_fam[['contig', 'bp']] = CSD_fam['merged'].apply(SNP_splitter)
    CSD_fam['family'] = fam
    CSD_fam.drop('merged', axis=1, inplace=True)
    # Split IDs to extract contig and base pair infos
    CSD_full = CSD_full.append(CSD_fam)

CSD_full.to_csv('data/SNP_lists/' + path.basename(thresh) +
'_raw_CSD_candidates.txt', index=False)


####################
# Visualizing data #
####################
for fam in sorted(fam_geno):
    fe = SNP_sum[fam]["Female het."]; me = SNP_sum[fam]["Male het."]
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
    if grouped == 'F':
        m1n = len(traits.Name[(traits.Ploidy=='H') & (traits.Family==fam)])
    else:
        m1n = len(traits.Name[traits.Ploidy=='H'])
        
    textstr = 'Threshold={0}\nM1N={1}\nM2N={2}\nF={3}\nSNPs={4}'.format(
        thresh.split('/')[-1], m1n, males[fam].shape[0], females[fam].shape[0],
        SNP_sum[fam].shape[0])
    plt.text(1.3, 1.2, textstr, fontsize=9)
    plt.draw()
    # plt.show()
    plt.savefig("reports/lab_book/assoc_explo_fam/" + fam + ".pdf")
