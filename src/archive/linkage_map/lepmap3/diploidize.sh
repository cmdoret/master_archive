# This script takes the STACKS populations output vcf file and a matching list of samples 
# with ploidy, sex and family informations as inputs. It outputs diploidized versions of 
# both files. 

# Here, diploidization is defined as:
# 1. Duplication of alleles for all haploid samples.
# 2. Duplications of all diploid samples.
# 3. Generation of a "synthetic father" by copying asexual mothers.

# Cyril Matthey-Doret
# 23.05.2018

### HANDLE ARGUMENTS ###

# Help message
function usage() {
    cat <<EOF
Usage: `basename $0` -v file.vcf -s samples.tsv [-h]
    -v  path to VCF file to diploidize.
    -s  path to tab-separated file containing the following columns with headers:
            Name: Alphanumeric value.
            Sex: Male or female [M,F].
            Family: Alphanumeric value.
            Generation: Parents or offspring [F3, F4]
            Ploidy: Haploid or diploid [H, D]
    -o  output folder where the updated vcf and sample list should be saved
    -h  display this help.
EOF
    exit 0
}

# Parsing CL arguments
while getopts ":v:s:o:h" opt
do
    case $opt in
        v ) VCF=${OPTARG} ;;
        s ) SPL=${OPTARG} ;;
        o ) OUT=${OPTARG} ;;
        h ) usage ;;
        \?) usage ;;
    esac
done

if [[ "x" == "x$VCF" ]] || [[ "x" == "x$SPL" ]] || [[ "x" == "x$OUT" ]]
then
    echo "Error: Input and output path must be provided."
    usage
    exit 0
fi


### SETUP ###

# Using this for optional loading bar. Ignore if absent.
source src/misc/jobs_manager.sh 2> /dev/null || true

# find 0-indexed indices of relevant variable to make script work even if 
# columns are swapped/added/removed
read -r -a headers < $SPL
for i in ${!headers[@]}
do
    case "${headers[$i]}" in
        'Name') NAME=$i
        ;;
        'Generation') GEN=$i
        ;;
        'Sex') SEX=$i
        ;;
        'Ploidy') PLOID=$i
        ;;
    esac
done

# Copying files to keep inputs unchanged
mkdir -p "$OUT"

OVCF="${OUT}/diploidized.vcf" 
cp "$VCF" "$OVCF"
OSPL="${OUT}/diploidized.tsv"
cp "$SPL" "$OSPL"

function get_sample() {
    # Given a sample ID, returns an int corresponding
    # to the (1-based) index of the column containing 
    # the sample genotype in the VCF file.

    local sample=$1
    # Line starting with a single hash in VCF (header line)
    header=($(grep -m 1 '^#[^#]' "$VCF"))
    
    # Loops over indices of columns (0-indexed)
    for i in "${!header[@]}";do
        if [[ "${header[$i]}" == "${sample}" ]]
        then
            # Return 1-indexed column number
            echo $((i+1))
        fi
    done
}

function allele_dup() {
    # Takes a column number and sample id as input.
    # Diploidize all genotypes in the vcf file at that 
    # column and edit the ploidy field for the corresponding 
    # sample in the sample list.

    local vcf_col=$1
    local spl_name=$2
    # At each position, duplicate allele of sample at input column
    awk -v id=$vcf_col 'BEGIN {OFS="\t"}
                        /^#/ {print}
                        /^[^#]/ {$id=$id"/"$id; print}' $OVCF > tmp_vcf.vcf \
            && mv tmp_vcf.vcf $OVCF

    # Edit ploidy in sample list
    awk -v name="$spl_name" \
        -v nameid=$((NAME+1)) \
        -v ploidid=$PLOID 'BEGIN {OFS="\t"}
                           $nameid == name {$ploidid="H"}{print}' $OSPL \
            > tmp_spl.tsv && mv tmp_spl.tsv "$OSPL"
}

function sample_dup() {
    # Takes an input column number and sample name and sex. 
    # Duplicates sample at target column. The new sample (of desired sex) #
    # is appended as the last column. The new sample is also appended to the sample list file

    local vcf_col=$1
    local spl_name=$2
    local spl_sex=$3
    
    # Appending new sample column in vcf file with modified name
    awk -v sname=$spl_name \
        -v ncol=$vcf_col 'BEGIN {OFS="\t"}
                          /^##/ {print $0}
                          /^#[^#]/ {print $0,sname"_DUP"}
                          /^[^#]/ {print $0,$ncol}' $OVCF > tmp_vcf.vcf \
        && mv tmp_vcf.vcf $OVCF

    # Appending new sample row in sample list with modified name and sex
    cat $OSPL <(grep -m 1 "$spl_name" "$OSPL" | \
        awk -v sexid=$((SEX+1))   -v sex=$spl_sex \
            -v nameid=$((NAME+1)) -v name=$spl_name 'BEGIN {OFS="\t"}
                {$sexid=sex;$nameid=name"_DUP";print}') \
            > tmp_spl.tsv && mv tmp_spl.tsv "$OSPL"
}

# Keep track of progress
tot=$(wc -l $SPL | awk '{print $1}')
curr=0

### RUN ###

# Loop over sample names
while read -a sample
do  
    # Loading bar or fallback echo
    prettyload $curr $tot 2> /dev/null || echo -ne "  $curr | $tot \r"
    ((curr++))
    col_ID=$(get_sample "${sample[$NAME]}")
    if [[ "${sample[$SEX]}" == 'F' ]]
    then
        # Mother: generate identical father
        if [[ "${sample[$GEN]}" == 'F3' ]]
        then
            sample_dup "$col_ID" "${sample[$NAME]}" "M"
        # Daughter: duplicate sample
        else
            sample_dup "$col_ID" "${sample[$NAME]}" "F"
        fi
    else
        # Haploid son: make it diploid
        if [[ "${sample[$PLOID]}" == 'H' ]]
        then
            allele_dup "$col_ID" "${sample[$NAME]}"
        # Diploid son: duplicate sample
        else
            sample_dup "$col_ID" "${sample[$NAME]}" "M"
        fi
    fi

done < <(tail -n +2 $SPL)
# Note we skip the header row when looping over samples
