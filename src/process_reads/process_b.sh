# This script processed the raw reads for a b-type library (adapter iA_12), using the process_radtags program from STACKS.
# Takes the library ID and optionally, the number of authorized mismatches in adaptor as an argument
# Cyril Matthey-Doret
# 15.03.2017

if [[ $# -eq 0 ]] ; then
    echo 'I need the ID of the library !'
    exit 0
fi
ID=$1
if [ $# -ge 2 ]
then
    MM=$1
else
    MM=2
fi

#BSUB -J radtags_b
#BSUB -o demulti-b-output_.o
#BSUB -e demulti-b-error_.e
##BSUB -n 4
#BSUB -u cmatthey@unil.ch
##BSUB -R "span[ptile=4]"
#BSUB -q normal
#BSUB -R "rusage[mem=4000]"
#BSUB -M 8000000

wd=/scratch/beegfs/monthly/cmatthey/data/processed
mkdir -p $wd'_'$MM'/lib'$ID

module add UHTS/Analysis/stacks/1.30;

process_radtags -p /scratch/beegfs/monthly/cmatthey/data/raw_reads/lib$ID/ \
-b /scratch/beegfs/monthly/cmatthey/data/barcodes/barcodes_radwasp$ID \
-o /scratch/beegfs/monthly/cmatthey/data/processed_$MM/lib$ID \
-r -q -c -e ecoRI --filter_illumina -i gzfastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm $MM;

module rm UHTS/Analysis/stacks/1.30;
