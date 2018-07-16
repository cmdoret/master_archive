#!/bin/bash
# Input: non-anchored (old) and anchored (new) genome assemblies, list of marker coordinates (in old assembly) and cM values
# Output: cM/BP ratio along genome, computed in all intervals between linkage map markers.
# Computes cM/bp relationship along genome using markers from Jens linkage map.
# Note: Marker must be sorted by ascending Chr,cM
# CMD, 20181607


usage(){

    cat <<USAGE
`basename $0` [OPTION]...
Compute cM/BP ratio between all markers in a linkage map and convert coordinates to an anchored assembly.
    -n, --new       path to the new (anchored) assembly fasta file.
    -r, --old_ref   path to the old (non-anchored) assembly fasta file.
    -m, --markers   list of markers, bp and cM values. 1 marker per line in the format contig,BP,cM.
    -p, --preserve  if available, use intermerdiate files from previous run.
    -o, --output    path to output file.
USAGE
    exit 0
}

while [[ $# -gt 0 ]]
do
key="$1"
case $key in
    -n|--new_ref)	NREF="$2";	shift;shift;;
    -r|--old_ref)	OREF="$2";	shift;shift;;
    -m|--markers)	MARK="$2";	shift;shift;;
    -o|--output)	OUTF="$2";	shift;shift;;
    -p|--preserve)  PRES="T";	shift;shift;;
    *)				usage;		shift;shift;; # Unknown option
esac
done

# Check if any argument is missing
if [ -z "$NREF" ];  then miss="New assembly"; fi
if [ -z "$OREF" ];  then miss="Old assembly"; fi
if [ -z "$MARK" ];  then miss="New assembly"; fi
if [ -z "$OUTF" ];  then miss="Output path";  fi

# Signal missing argument, print help and exit
if [ ! -z ${miss+x} ];then echo "$miss not provided."; usage; fi

# Compute coordinates of markers in new assembly and store correspondance list to file
# Unless preserve was specified (in which case, an existing file will be correspondance list will be used)

if [ ! -z ${PRES+x} ] || [ ! -f "$OUTF.corresp" ]; then
    rm "$OUTF.corresp"
    while read marker; do
        # For all markers: looking for surrounding sequence in newer assembly 
        python2 src/convert_coord/contig2chr.py "$OREF" "$NREF" --pos1 "$marker" \
            | cut -f1,2 -d ',' \
            | paste -d ',' <(echo $marker) - \
            | grep 'chr' >> "$OUTF.corresp"
done < "$MARK"
fi

# Compute average cM/BP per chromosome
awk 'BEGIN{FS=",";OFS=","}
     { $4 == chr {BP=$4; if ($3 > cMax) cMax=$3} 
       $4 != chr {print chr,(cMax/BP); chr=$4; cMax=$3}
     }' "$OUTF.corresp" > "$OUTF.cMean"

prev=(0 0 0 0 0)
rm "$OUTF.csv"
while read -a m; do
    # If still on same contig, compute bp/cM ratio between previous and current marker
    if [ ${prev[0]} == ${marker[0]} ] || [ ${prev[0]} == 0 ]; then
        cM=$(echo "${marker[2]} - ${prev[2]}" | bc -l)
        bp=$(( ${marker[4]} - ${prev[4]} ))
        ratio=$(echo "$cM / $bp"  | bc -l | xargs printf "%.3f\n")
    else
        # if still on same chromosome, use mean chrom. ratio to estimate inter-contig ratio
        if [ ${prev[3]} == ${prev[3]} ]; then
            ratio=$(awk -v chrom="${marker[3]}" '{$1 == chrom {print $2}}' "$OUTF.cMean")
        else
            # New chromosome. No need to subtract (previous is 0)
            ratio=$(echo "${marker[2]} / ${marker[4]}" | bc -l | xargs printf "%.3f\n")
            prev[4]=0
        fi
    fi
    # send chromosome, interval between previous and current markers in BP, and cM/BP ratio in interval
    echo "${marker[3]},${prev[4]},${marker[4]},$ratio" >> "$OUTF.csv"
    prev=(${marker[@]})
done < "$OUTF.corresp"
