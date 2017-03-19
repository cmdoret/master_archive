#BSUB -J post_qc
#BSUB -o post_qc.o
#BSUB -e post_qc.e
#BSUB -q normal

lib=7
MM=2

module add UHTS/Quality_control/fastqc/0.11.2
fastqc -o ./data/processed_$MM/qc/ ./data/processed_$MM/lib$lib/lib$lib.fq
module rm UHTS/Quality_control/fastqc/0.11.2 
