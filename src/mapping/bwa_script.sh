

data_dir=data/processed/
index=data/ref_genome/ordered_genome/merged.fasta ## path and prefix of indexed genome files

threads=2

# Parameters I want to change
ALG=mem ## mem or backtrack[mem/aln]
MM=4 ## mismatches (aln)[4]
K=19 ## min seed length (mem)[19]
W=100 ## band width (mem)[100]
prefix=BWA
out_dir=data/mapped/$ALG-$MM-$K-$W/

## create output directory for bam files:
mkdir -p $out_dir/bam

for sample in $(cut -f1 data/popmap.tsv) #this is the list of sample names
do
  # Limit number of queued jobs to 100 at a time
  while [ $(bjobs -w | grep 'bam' | wc -l) -gt 100 ]
  do
    sleep 1;
  done
  echo "processing sample $sample";
  # Sending each sample as a separate jobs
  bsub <<MAPSAMPLE
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
    
    # select aln or mem
    case $ALG in
        'mem')
            bwa mem -k $K -t $threads -w $W  $index $data_dir/$sample.fq.gz > $out_dir/$sample-$prefix.sam
            ;;
        'aln')
            bwa aln -n $MM -t $threads $index $data_dir/$sample.fq.gz > $out_dir/$sample.sai ## align reads
            bwa samse -n 3 $index $out_dir/$sample.sai $data_dir/$sample.fq.gz > $out_dir/$sample-$prefix.sam # index alignment file
            ;;
    esac

    # perl script removes reads which map more than once
  	perl src/mapping/split_sam.pl -i $out_dir/$sample-$prefix.sam -o $out_dir/$sample-$prefix >> $out_dir/split_summary.log
    
    # Remove original SAM files
    rm -v $out_dir/$sample-$prefix.sam
    
    # Convert SAM files to BAM
  	samtools view -@ $threads -bS -o $out_dir/bam/$sample-$prefix-uniq.bam $out_dir/$sample-$prefix-uniq.sam
    
    # Sort alignments by leftmost coordinate
    samtools sort -@ $threads $out_dir/bam/$sample-$prefix-uniq.bam -o $out_dir/bam/$sample-$prefix-uniq.sorted.bam
    
    # Index BAM files
    samtools index $out_dir/bam/$sample-$prefix-uniq.sorted.bam
    
    # Output index statistics
  	samtools idxstats $out_dir/bam/$sample-$prefix-uniq.sorted.bam
    
    # Remove unsorted bam files
  	rm -v $out_dir/bam/$sample-$prefix-uniq.bam
  	date
    
    # Compress SAM files
    gzip -v $out_dir/$sample*.sam
MAPSAMPLE
done

# Wait for all mapping jobs to be finished before resuming pipeline
while [ $(bjobs -w | awk '/RUN/ {print $7}' | grep 'bam' | wc -l) -gt 0 ]
do
    sleep 2;
done
