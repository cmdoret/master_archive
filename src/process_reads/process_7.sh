# This script processed the raw reads in the lib7 folder, using the process_radtags program from STACKS.
# Takes the number of authorized mismatches in adaptor as an argument
# Cyril Matthey-Doret
# 15.03.2017

if [ $# -ge 1 ]
then
    MM=$1
else
    MM=2
fi

#BSUB -J radtags7
#BSUB -o demulti-7-output.o
#BSUB -e demulti-7-error.e
##BSUB -n 4
#BSUB -u cmatthey@unil.ch
##BSUB -R "span[ptile=4]"
#BSUB -q normal
#BSUB -R "rusage[mem=4000]"
#BSUB -M 8000000

wd=/scratch/beegfs/monthly/cmatthey/data/processed
mkdir -p $wd'_'$MM'/lib7'

module add UHTS/Analysis/stacks/1.30;

process_radtags -p /scratch/beegfs/monthly/cmatthey/data/raw_reads/lib7/ \
-b /scratch/beegfs/monthly/cmatthey/data/barcodes/barcodes_radwasp7 \
-o /scratch/beegfs/monthly/cmatthey/data/processed_$MM/lib7 \
-r -q -c -e ecoRI  --filter_illumina -i gzfastq --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG --adapter_mm $MM;

module rm UHTS/Analysis/stacks/1.30;
