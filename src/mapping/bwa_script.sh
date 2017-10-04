

#module add UHTS/Aligner/bowtie2/2.2.4;
data_dir=data/processed/
index=data/ref_genome/ordered_genome/merged.fasta ## path and prefix of indexed genome files

threads=12

# Parameters I want to change
ALG=mem ## mem or backtrack[mem/aln]
MM=4 ## mismatches (aln)[4]
K=19 ## min seed length (mem)[19]
W=100 ## band width (mem)[100]
#D=100 ## max dist between query and ref position (mem)
#R=1.5 ## min reseeding length (mem)
prefix=BWA
out_dir=data/mapped/$ALG-$MM-$K-$W/

## Do some work:
mkdir -p $out_dir/bam

date

for sample in $(cut -f1 data/popmap.tsv) #this is a list of sample names (similar as the popmap file for Populations, but without reproductive mode)
do
  while [ $(bjobs -w | awk '/RUN/ {print $7}' | grep 'bam' | wc -l) -gt 40 ]
  do
    sleep 1;
  done
  echo "\nprocessing sample $sample\n";
  #cp -v $data_dir/$sample* $sample.fq.gz
  #gunzip -v $sample.fq.gz
  bsub <<TMPSCRIPT
    #!/bin/bash

    #BSUB -L /bin/bash
    #BSUB -o BWA-%J-OUT.txt
    #BSUB -e BWA-%J-ERROR.txt
    #BSUB -u cmatthey@unil.ch
    #BSUB -J bam
    #BSUB -n 2
    #BSUB -R "span[ptile=2]"
    #BSUB -q priority
    #BSUB -R "rusage[mem=4000]"
    #BSUB -M 4000000

    module add UHTS/Aligner/bwa/0.7.2;
    module add UHTS/Analysis/samtools/1.3;

    case $ALG in
        'mem')
            bwa mem -k $K -t $threads -w $W  $index $data_dir/$sample.fq.gz > $out_dir/$sample-$prefix.sam
            ;;
        'aln')
            bwa aln -n $MM -t $threads $index $data_dir/$sample.fq.gz > $out_dir/$sample.sai ## align reads
            bwa samse -n 3 $index $out_dir/$sample.sai $data_dir/$sample.fq.gz > $out_dir/$sample-$prefix.sam # index alignment file
            ;;
    esac

  	perl src/mapping/split_sam.pl -i $out_dir/$sample-$prefix.sam -o $out_dir/$sample-$prefix >> $out_dir/split_summary.log ## change perl script path (script removes reads which map more than once)
  	rm -v $out_dir/$sample-$prefix.sam
      # cd $out_dir/bam/
  	samtools view -@ $threads -bS -o $out_dir/bam/$sample-$prefix-uniq.bam $out_dir/$sample-$prefix-uniq.sam ## converts from SAM to BAM
  	samtools sort -@ $threads $out_dir/bam/$sample-$prefix-uniq.bam -o $out_dir/bam/$sample-$prefix-uniq.sorted.bam ## sort bam
  	samtools index $out_dir/bam/$sample-$prefix-uniq.sorted.bam # index bam
  	samtools idxstats $out_dir/bam/$sample-$prefix-uniq.sorted.bam # get index stats
  	rm -v $out_dir/bam/$sample-$prefix-uniq.bam ## remove unsorted bam file
  	date
    gzip -v $out_dir/$sample*.sam # gzip sam files
TMPSCRIPT
done



while [ $(bjobs -w | awk '/RUN/ {print $7}' | grep 'bam' | wc -l) -gt 0 ]
do
    sleep 2;
done
