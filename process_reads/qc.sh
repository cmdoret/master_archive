 
#BSUB -J qc
#BSUB -o qc.o
#BSUB -e qc.e
#BSUB -q normal

lib=7

module add UHTS/Quality_control/fastqc/0.11.2
fastqc -o ./data/raw_reads/qc/ ./data/raw_reads/lib$lib/lib$lib.fq
module rm UHTS/Quality_control/fastqc/0.11.2
