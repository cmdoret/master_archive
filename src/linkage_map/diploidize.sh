# This script takes the STACKS populations output vcf file and a matching list of samples 
# with ploidy, sex and family informations as inputs. It outputs diploidized versions of 
# both files. 
# Here, diploidization is defined as:
# 1. Duplication of alleles for all haploid samples.
# 2. Duplications of all diploid samples.
# 3. Generation of a "synthetic father" by copying asexual mothers.
# Cyril Matthey-Doret
# 23.05.2018

# Help message
function usage() {
    cat <<EOF
Usage: `basename $0` ...
    -...    ...
EOF
    exit 0
}

# Parsing CL arguments
while getopts ":v:s:h" opt
do
    case $opt in
        v ) VCF=${OPTARG} ;;
        s ) SPL=${OPTARG} ;;
        h ) usage ;;
        \?) usage ;;
    esac
done

function get_sample() {
    # Given a vcf file and a sample ID, returns an int corresponding
    # column containing the sample genotype
    local vcf=$1
    local sample=$2
    # Line starting with a single hash in VCF (header line)
    header=($(grep -c 1 '^#[^#]' $vcf))
    
    # Loops over indices of columns
    for i in "${!header[@]}";do
        if [[ "${my_array[$i]}" = "${sample}" ]]:then
            return "${i}";fi
    done
}

function allele_dup() {

    
}

function sample_dup() {
}

# Loop over sample names
while read sample
do
    col_ID=$(get_sample $sample)
    if sex == 'F'
        # If mother: generate identical father
        if Generation == 'F3'
            sample_dup M $col_ID "${sample}_DUP"
        # If daughter duplicate sample
        else
            sample_dup F $col_ID "${sample}_DUP"
    else
        # If haploid son: make it diploid
        if Ploidy == 'H'
            allele_dup $col_ID
        # If diploid son: duplicate sample
        else
            sample_dup M $col_ID "${sample}_DUP"

done < $SPL
