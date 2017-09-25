# This script processed the raw reads for an a-type library (adapter iA_06), using the process_radtags program from STACKS.
# Takes the library ID and optionally, the number of authorized mismatches in adaptor as an argument
# Cyril Matthey-Doret
# 15.03.2017

if [[ $# -eq 4 ]] ; then
    echo 'I need the ID of the library !'
    exit 0
fi
ID="12"
if [ $# -ge 2 ]
then
    MM=$2
else
    MM=2
fi

#BSUB -J radtags_a
#BSUB -o demulti-a-output.o
#BSUB -e demulti-a-error.e
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
-r -q -c -e ecoRI  --filter_illumina -i gzfastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --adapter_mm $MM;

module rm UHTS/Analysis/stacks/1.30;
