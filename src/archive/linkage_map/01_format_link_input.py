# Centralizes all of the formatting required before running LEPMAP3
# Given the following files:
# 1. List of sequenced samples (fixed.tsv file)
# 2. List of mother-offspring correspondance (families.tsv file)
# 3. A VCF file of all sequenced samples
# It will:
# 1. Reconstruct the genotypes of mothers from their offspring
# 2. duplicate the mothers to generate artificial fathers
# 3. Merge pairs of haploid offspring into synthetic diploid and remove unpaired haploids
# 4. Format all of these processed individuals into a LEPMAP compatible pedigree and a VCF file.

# Using scikit-allel to parse VCF
import  allel
import numpy as np
import pandas as pd

#### TODO: Parse CLI arguments + help function

#### DEFINITIONS


def append_mother(df):
	"""
	Appends a synthetic mother entry to the dataframe if no mother is present.

    Args:
        df: pandas dataframe containing members of a family.

    Returns:
        pandas.DataFrame : Input dataframe with an additional entry (row)
	"""
	# If no mother present in the family: add one
	if not (df.Generation=="F3").any():
		nsamples = df.shape[0]
		df = df.append(df.iloc[1,:], ignore_index=True)
		df.loc[nsamples,:] = "NA"
		df.loc[nsamples,["Name", "Sex", "Family", "Generation", "Ploidy"]] = \
		df.loc[0,"Parent_id"], "F", df.loc[0,"Family"], "F3", "D"
	return df



def reconstruct_mother():
	pass


def find_pairs():
	pass

#### READ

# Use scikit-allel to parse vcf file and build genotype array
in_vcf = allel.read_vcf("data/linkage_map/lib13/grouped_d-5_r-80/batch_1.vcf")
in_gt = allel.GenotypeArray(in_vcf['calldata/GT'])

# Loade ploidy and family files
ploid = pd.read_csv("data/linkage_map/lib13/fixed.tsv", sep="\t")
fam = pd.read_csv("data/families.tsv", sep="\t", )

#### FORMAT INPUT

fam = fam.loc[:,["Number", "Parent_id"]]
in_vcf["samples"]
# Combine both sample description tables into a single one
indv = ploid.merge(fam, left_on='Name', right_on='Number', \
				   how='left').drop('Number', axis=1)
# Append artificial mothers at the end of each family
indv = indv.groupby("Family").apply(append_mother)
# Remove grouping factor
indv.reset_index(drop=True, inplace=True)

# Sample name to index correspondance in VCF file
sample_dict ={v:i for i,v in enumerate(in_vcf["samples"])}
vcf_corr_id = pd.DataFrame.from_dict(sample_dict, orient='index', \
									 columns=["vcf_idx"])
indv = indv.merge(vcf_corr_id, left_on='Name', right_index=True, how='left')

# Making subpops by family in vcf
vcf_families = {}
for fam in indv.Family.unique():
	fam_idx = indv.loc[(indv.Family == fam) & (indv.Generation == "F4"),"vcf_idx"]
	vcf_families[fam] = fam_idx.astype('int32').tolist()

#### RECONSTRUCT Mothers
# Compute allele frequencies in each family
fam_freq = in_gt.count_alleles_subpops(vcf_families)
# Compute allele frequencies in family at each position

list(map(lambda x: fam_freq[x]/fam_freq[x].sum(1)[:,None], fam_freq))

#### GENERATE FATHERS

#### MERGE HAPLOIDS

#### FORMAT OUTPUT
