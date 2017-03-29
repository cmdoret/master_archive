#!/bin/bash

data_dir=../../data/processed/
index=../../data/ref_genome/lfabarum_bwa_index ## path and prefix of indexed genome files

threads=4

# Parameters I want to change
ALG=mem ## mem or backtrack[mem/aln]
MM=4 ## mismatches (aln)[4]
K=19 ## min seed length (mem)[19]
W=100 ## band width (mem)[100]
#D=100 ## max dist between query and ref position (mem)
#R=1.5 ## min reseeding length (mem)
prefix=BWA
out_dir=/scratch/beegfs/monthly/cmatthey/data/mapped/$ALG-$MM-$K-$W/

## Do some work:
mkdir -p bam

date

for sample in $(tail -n +2 /scratch/beegfs/monthly/cmatthey/data/F4_table.csv | awk -F, '/01/ {print $1"_"$2}') #this is a list of sample names (similar as the popmap file for Populations, but without reproductive mode)
do 
        echo "\nprocessing sample $sample\n";
        #cp -v $data_dir/$sample* $sample.fq.gz
        #gunzip -v $sample.fq.gz
        
        case $ALG in
            'mem')
                bwa mem -k $K -t $threads -w $W  $index $data_dir/$sample.fq.gz > $sample.sai
                ;;
            'aln')
                bwa aln -n $MM -t $threads $index $data_dir/$sample.fq.gz > $sample.sai ## align reads
                ;;
        esac
        
        bwa samse -n $MM $index $sample.sai $data_dir/$sample.fq.gz > $sample-$prefix.sam # index alignment file
	perl /scratch/beegfs/monthly/cmatthey/src/mapping/split_sam.pl -i $sample-$prefix.sam -o $sample-$prefix >> split_summary.log ## change perl script path (script removes reads which map more than once)
	rm -v $sample-$prefix.sam
	cd bam/
	samtools view -@ $threads -bS -o $sample-$prefix-uniq.bam ../$sample-$prefix-uniq.sam ## converts from SAM to BAM
	samtools sort -@ $threads $sample-$prefix-uniq.bam -o $sample-$prefix-uniq.sorted.bam ## sort bam
	samtools index $sample-$prefix-uniq.sorted.bam # index bam
	samtools idxstats $sample-$prefix-uniq.sorted.bam # get index stats
	rm -v $sample-$prefix-uniq.bam ## remove unsorted bam file
	cd ../
	gzip -v *.sam # gzip sam files
	mv *.sam.gz bam/*.bam $out_dir
	date
done
