# This script processed the raw reads in the lib7 folder, using the process_radtags program from STACKS.
# Takes number of mismatches in adaptor as an argument.
# Cyril Matthey-Doret
# 15.03.2017
#BSUB -J radtags7b
#BSUB -o 'demulti-%J-output_'$MM'.txt'
#BSUB -e 'demulti-%J-error_'$MM'.txt'
#BSUB -n 12
#BSUB -u cmatthey@unil.ch
#BSUB -R "span[ptile=4]"
#BSUB -q long
#BSUB -R "rusage[mem=4000]"
#BSUB -M 8000000

if [ $# -ge 1 ]
then
    MM=$1
else
    MM=2
fi

wd=/scratch/beegfs/monthly/cmatthey/data/processed
mkdir -p $wd'_'$MM'/lib7b'

module add UHTS/Analysis/stacks/1.30;

process_radtags -p /scratch/beegfs/monthly/cmatthey/data/raw_reads/lib7b/ \
-b /scratch/beegfs/monthly/cmatthey/data/barcodes/barcodes_radwasp7b \
-o /scratch/beegfs/monthly/cmatthey/data/processed_$MM/lib7b \
-r -q -c -e ecoRI --filter_illumina --adapter_1 GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG --adapter_mm $MM;

module rm UHTS/Analysis/stacks/1.30;
