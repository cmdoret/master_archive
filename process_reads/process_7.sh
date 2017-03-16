# This script processed the raw reads in the lib7 folder, using the process_radtags program from STACKS.
# Cyril Matthey-Doret
# 15.03.2017

#BSUB -o demulti-%J-output.txt
#BSUB -e demulti-%J-error.txt
#BSUB -J radtags7
#BSUB -n 12
###BSUB -u cmatthey@unil.ch
###BSUB -R "span[ptile=4]"
##BSUB -q long
#BSUB -R "rusage[mem=4000]"
#BSUB -M 8000000


module add UHTS/Analysis/stacks/1.30;

process_radtags -p /scratch/beegfs/monthly/cmatthey/data/raw_reads/lib7/ \
-b /scratch/beegfs/monthly/cmatthey/data/barcodes/barcodes_radwasp7 \
-o /scratch/beegfs/monthly/cmatthey/data/processed/lib7 \
-r -q -c -e ecoRI --filter_illumina --adapter_1 GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT --adapter_mm 2;

module rm UHTS/Analysis/stacks/1.30;
