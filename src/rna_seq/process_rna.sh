#!/bin/bash
# This script takes a sorted bam file containing reference-aligned RNA-seq reads
# And the the reference genome in FASTA format as input.
# Cyril Matthey-Doret
# 11.11.2017
#BSUB -J RNA_assembly
#BSUB -q normal
#BSUB -n 36
#BSUB -M 32000000
#BSUB -R "span[ptile=36]"
#BSUB -R "rusage[mem=32000]"

module add UHTS/Assembler/cufflinks/2.2.1
module add UHTS/Analysis/samtools/1.3

rna='../../data/rna_seq/'
ref='../../data/ref_genome/canu2_low_corrRE.contigs.fasta'
align="${rna}STAR_larvae_aligned.sort.bam.bam"

rm -rf "${rna}assembled/"
mkdir -p "${rna}assembled"

# Assemble transcriptome using cufflinks
cufflinks -p 32 -b $ref -u --library-type ff-firststrand  $align -o "${rna}assembled/"

# Generate bed file with transcripts positions
