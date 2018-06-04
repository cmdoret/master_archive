# This script takes a VCF file  and a list of samples. It uses vcftools to filter 
# samples that will be used to build the linkage map. It then runs lepmap3 to build 
# a linkage map.
# Cyril Matthey-Doret
# 24.05.2018

# Help message
function usage() {
    cat <<EOF
Usage: `basename $0` -v file.vcf -p pedigree.tsv [-h]
    -v  input VCF file (must have been diploidized using (diploidize.sh).
    -p  pedigree file containing samples that must be used for the linkage map.
    -o  output directory.
    -h  display this help.
EOF
    exit 0
}

# Parsing CL arguments
while getopts ":v:p:o:h" opt
do
    case $opt in
        v ) VCF=${OPTARG} ;;
        p ) PED=${OPTARG} ;;
        o ) OUT=${OPTARG} ;;
        h ) usage ;;
        \?) usage ;;
    esac
done

if [[ "x" == "x$VCF" ]] || [[ "x" == "x$PED" ]] || [[ "x" == "x$OUT" ]]
then
    echo "Error: VCF, Pedigree and output path must be provided."
    usage
    exit 0
fi

# Copying files to keep inputs unchanged
mkdir -p "$OUT"
VCF2="${OUT}/linkage_samples.vcf"
source src/misc/dependencies.sh

#  Filter samples in VCF file using those in the pedigree file
samples=$(cut -f3- "$PED" | sed -n '2p')

vcftools --keep <(echo "$samples" | tr '[\t ]' '\n') \
         --vcf "$VCF" \
         --recode \
         --out "$VCF2"

mv "$VCF2.recode.vcf" "$VCF2"

# Calling genotypes and removing uninformative markers
java -classpath "$LEPMAP3" ParentCall2 data="$PED" \
                                       vcfFile="$VCF2" \
                                       removeNonInformative=1 > "$OUT/01_geno.call"

# Filtering out markers  with too many missing calls or high segregation disortion
java -classpath "$LEPMAP3" Filtering2 data="$OUT/01_geno.call" \
                                      dataTolerance=0.001 > "$OUT/02_filtered.call"

# Assigns markers into linkage groups (LGs)
java -classpath "$LEPMAP3" SeparateChromosomes2 data="$OUT/02_filtered.call" \
                                                lodLimit=2 \
                                                sizeLimit=10 \
                                                numThreads=4 \
                                                distortionLod=1 > "$OUT/03_linkmap.txt"

# Tries to join additional single markers to LGs, by computing LOD score between singles and LGs
java -classpath "$LEPMAP3" JoinSingles2All map="$OUT/03_linkmap.txt" \
                                           data="$OUT/02_filtered.call" \
                                           numThreads=4 \
                                           lodLimit=2 > "$OUT/04_linkmap_js.txt"

# Find number of LGs
N_LG=$(cut -f1 "$OUT/04_linkmap_js.txt" | sort -n | uniq | tail -n 1)

# Order markers on each LG to optimise likelihood 
for (( LG=1; LG<=$N_LG; LG++ ))
do
    java -classpath "$LEPMAP3" OrderMarkers2 map="$OUT/04_linkmap_js.txt" \
                                             data="$OUT/02_filtered.call" \
                                             randomPhase=1 \
                                             numThreads=4 \
                                             chromosome=$LG > "$OUT/05_order_LG$LG.txt"
done
