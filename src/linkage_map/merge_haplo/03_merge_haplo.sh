#!/usr/bin/env bash
# This script is used to merge haploid into diploid samples for LEPMAP3
# It operates on ParentCall2 output (because it is easier than on VCF files).
# CMD
# 18.06.2018

# In file, genotypes are in following order:
# 		AA AC AG AT CC CG CT GG GT TT
# Use index of maximum value to get haplotype
# 1:A 2:C 4:G 8:T 0:N
# Note N should also be for cases where there is no single max value
# Encode each possible haplotype as binary for easy sum afterwards
#      N A N N N C N N G N T
g2bin=(0 1 0 0 0 2 0 0 4 0 8)

# Generate pairs of samples among offspring
python3 define_pairs.py

# Stores 1-indexed column number of samples
# Pair # N consists of ${to_keep[$N]} and ${to_drop[$N]}
# In each pair, the haplotype of to_drop will be used to make
# to_keep diploid. At the end, all to_drop are removed
to_keep=($(sed -n '1p' sample_pairs.txt))
to_drop=($(sed -n '2p' sample_pairs.txt))
# Samples to exclude. Consists of parents and unpaired offpsring (if odd number in family)
to_excl=($(sed -n '3p' sample_pairs.txt))

function get_geno {
# Get genotype with highest likelihood for a given sample
		sample_idx="$1"
		genotypes="$2"
    # Extract sample genotype likelihoods
		# Iterate through likelihoods to find actual genotype
    geno=$(echo "$genotypes" \
		| awk -v idx_g="$sample_idx" \
					 'BEGIN{FS="\t";OFS=" "}
					 {
						 print $idx_g
					 }' \
		| awk 'BEGIN{FS=" ";iMAX=0}
					{
						for (i=1;i<=NF;i++)
						{
							if ($i > MAX) {MAX=$i;iMAX=i}
							else if ($i == MAX) SCD=$i
						}
					}
					END{
						if (MAX==SCD) print 0;
					else print iMAX
					}')
    # If no max likelihood return 0, else binary code of haplotype

		echo ${g2bin[$geno]}
}

function set_geno {
# Write new genotype by combining two haplotypes
# Args: geno1 geno2 keep_index outpath
	geno_k="$1"
	geno_d="$2"
	keep_idx="$3"
	# If at least 1 missing / wrong haplotype: missing
	if [[ "$geno_k" -eq 0 ]] || [[ "$geno_d" -eq 0 ]];then
		diplo_gen=0
	# If 2 correct haplotypes: homozygous/heterozygous
	else
		diplo_gen=$(($geno_k+$geno_d))
	fi
	# Make sample index 0 indexed for array access
	keep_idx=$((keep_idx-1))
	format_geno=(0 0 0 0 0 0 0 0 0 0)
	if [[ $diplo_gen -gt 0 ]];then
		format_geno[$((diplo_gen-1))]=1
	fi
	#echo "$keep_idx"
	#echo "${format_geno[@]}"
	snp_write[$keep_idx]=$(echo -n "${format_geno[@]}" | sed 's/ /,/g')
}

outfile='test.txt'
rm "tmp_outfile.txt"
i=0
while read snp
do
		# Storing a copy of the record in memory for editing.
		snp_write=($(echo "$snp" | sed 's/ /,/g'))
    for n in ${!to_keep[@]}
    do
				# Get binary encoding of haplotype
        bin_k=$(get_geno ${to_keep[$n]} "$snp")
        bin_d=$(get_geno ${to_drop[$n]} "$snp")

				# Set genotype of merged diploid from sum of haplotype binary codes
				set_geno $bin_k $bin_d "${to_keep[$n]}"
    done
		# Writing record (SNP) to temporary output file
		# Substituting spaces for tabs, and then commas for spaces to restore
		# original formatting
		echo "${snp_write[@]}" \
		| sed 's/ /	/g' \
		| sed 's/,/ /g' >> "tmp_outfile.txt"
	((i++))
done < <(tail -n +7 01_geno.call)

# Include pedigree header in output
head -n 8 "01_geno.call" > $outfile
# Include only genotypes of samples to keep
cut -d"\t" -f1,2,$(echo -n ${to_keep[@]} | sed 's/ /,/g') "tmp_outfile.txt" >> $outfile


# Set index geno1 + geno2 to 1, others to 0
# Delete second sample of each pair
