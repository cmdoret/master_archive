# This script centralizes dependency path for the whole project
# Cyril Matthey-Doret
# 02.02.2018


# Loaded from LFS system
module add UHTS/Aligner/bwa/0.7.15 \
           UHTS/Analysis/samtools/1.4 \
           UHTS/Analysis/vcftools/0.1.15 \
           R/3.4.2 \
           Blast/ncbi-blast/2.6.0+ \
           UHTS/Analysis/BEDTools/2.26.0 \
           UHTS/Quality_control/fastqc/0.11.5 \
           UHTS/Analysis/deepTools/2.5.4 \
           SequenceAnalysis/GenePrediction/augustus/3.2.3 \
           SequenceAnalysis/HMM-Profile/hmmer/3.1b2 \
           UHTS/Analysis/freebayes/1.0.0 \
           UHTS/Analysis/stacks/1.48 \
           UHTS/Analysis/htslib/1.6

# Loaded from local install
export PATH=$PATH:"~/scratch/softwares/MCScanX/"
export PATH=$PATH:"~/scratch/softwares/busco/3.0.2b/scripts/"
export PATH=$PATH:"~/scratch/softwares/cufflinks-2.2.1/"

# Necessary for VCFtools perl scripts
export PERL5LIB=/software/UHTS/Analysis/vcftools/0.1.15/lib64/perl5/site_perl/5.24.0/

# AUGUSTUS config required to use BUSCO
export AUGUSTUS_CONFIG_PATH=/path/to/your/augustus_config
