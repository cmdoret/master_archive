# Transforms basepairs into cM in the case_control_all.tsv file 
# from the association mapping using the output of bp2cm.sh
# Assumes both files are sorted by ascending Chr, BP
# TODO: distance between pairs of SNPs needs to be ratio*D
# Where D is the distance between a SNP and end of interval.
# If SNPs on different intervals, compute D per interval and sum them
# TODO: increment cumulative sum after every pair of SNP to get actual cM unit.

#Start at line 2 to skip header
n_interv=2
echo -e "$(head -n 1 $hits)\tcM" > $out_hits

# Loop over association mapping hits
while read -a hit; do
    # Loop over linkage intervals 
    # Stop reading if interval starts after hit or is on diff chrom.
    while read -a interv && [ ${interv[0]} == ${hit[1]} ] \
                         && [ ${interv[1]} -gt ${hit[2]} ]; do
        # If hit is inside segment (i.e. before its end)
        if [ ${hit[3]} -le ${inter[2]}]; then
            echo -e "${hit[*]}\t${inter[2]}" >> "$out_hits"
        fi
        # SNP not found, next read will not use this interval
        ((n_interv++))
    # Only check intervals from n_inter to the end
    done < <(sed -n "${n_inter},p" $bp2cm)
done < $hits
