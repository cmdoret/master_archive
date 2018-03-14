#!/bin/bash
# This script handles the mapping of whole genome
# sequencing illumina reads from wild mothers.
# Cyril Matthey-Doret
# 11.10.2017


# Parsing CL arguments
while [[ "$#" > 1 ]];
do
    case $1 in
        # Working directory. Must contain a "raw" folder with reads
        --workdir) wgs="$2";;
        # Reference genome. Must be indexed beforehand
        --ref) index="$2";;
        *) break;;
esac; shift; shift
done

# #Mismatches allowed per read for mapping
MM=4
# #threads used when mapping
threads=8
# Max number of parallel jobs
MAXPROC=30

# Allows to cap # parallel jobs and wait for their execution
source src/misc/jobs_manager.sh

in_dir="${wgs}/raw/"

# Each sample is split into forward and reverse on 2 different lanes
# Building array of sample names from raw reads files
samples=( $(find "$in_dir" | grep "fastq" | sed 's/.*\(HYI-[0-9]*\)_R.*/\1/g' | sort | uniq) )

# Reinitializing folders
rm -rf "${wgs}/merged/"
mkdir -p "${wgs}/merged/"

tmp_dir="${wgs}/tmp/"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"

map_dir="${wgs}/mapped/"
rm -rf "$map_dir"
mkdir -p "$map_dir"

logs="${wgs}/log"
rm -rf "$logs"
mkdir -p "$logs"

for sample in ${samples[@]};
do

# Hang script if too many parallel jobs running
bmonitor "WGSBWA" $MAXPROC

# Submit the mapping of each sample as an independent job
bsub << EOF
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o ${logs}/${sample}-OUT.txt
#BSUB -e ${logs}/${sample}-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSBWA-${sample}
#BSUB -n 8
#BSUB -R "span[ptile=8]"
#BSUB -q normal
#BSUB -R "rusage[mem=64000]"
#BSUB -M 64000000

# Loading softwares
source src/misc/dependencies.sh

forward="${wgs}/merged/${sample}_R1.fastq.gz"
reverse="${wgs}/merged/${sample}_R2.fastq.gz"

# Merge lanes, keep forward and reverse separated
find "${in_dir}" -name "*${sample}*R1*" -type f | \
  sort | xargs cat > "\$forward"

find "${in_dir}" -name "*${sample}*R2*" -type f | \
  sort | xargs cat > "\$reverse"

# Trimmed reads filenames

# Paired:
trimF=\$(echo \$forward | sed 's%/\([^/]*\)$%/trim.\1%')
trimR=\$(echo \$reverse | sed 's%/\([^/]*\)$%/trim.\1%')

# Unpaired:
trimFU=\$(echo \$trimF | sed 's/R1.fastq.gz$/R1U.fastq.gz/')
trimRU=\$(echo \$trimR | sed 's/R2.fastq.gz$/R2U.fastq.gz/')

# Trim low quality ends
trimmomatic PE \$forward \$reverse \$trimF \$trimFU \$trimR \$trimRU \
            -trimlog ${logs}/${sample}-trim.log \
            LEADING:20 TRAILING:20

# Mapping paired ends reads
bwa mem -M -t $threads $index \$trimF \$trimR > "${map_dir}/${sample}.sam"

# Convert SAM files to BAM
samtools view -@ $threads -b -o "${map_dir}/${sample}.bam" "${map_dir}/${sample}.sam"

# Sort alignments by read name
samtools sort -@ $threads -n "${map_dir}/${sample}.bam" -o "${map_dir}/${sample}.nsort.bam"

# Fix mate information (adds ms and MC tags for markdup)
samtools fixmate "${map_dir}/${sample}.nsort.bam" "${map_dir}/${sample}.fixed.bam"

# Sort alignments by read name
#samtools sort -@ $threads -n "${map_dir}/${sample}.fixed.bam" -o "${map_dir}/${sample}.fixed.nsort.bam"

# Remove PCR duplicates (input needs to be name-sorted, otherwise secondary
# alignments are not considered in the duplicates)
#picard MarkDuplicates \
#      I="${map_dir}/${sample}.nsort.bam" \
#      O="${map_dir}/${sample}.dedup.bam" \
#      M="${map_dir}/${sample}.dup_metrics.txt" \
#      ASSUME_SORT_ORDER=queryname \
#      REMOVE_DUPLICATES=true

# Sort alignments by leftmost coordinate
samtools sort -@ $threads "${map_dir}/${sample}.fixed.bam" \
              -o "${map_dir}/${sample}.fixed.csort.bam"
# Index BAM files
samtools index "${map_dir}/${sample}.fixed.csort.bam"

# Remove sam and unsorted/temporary bam files
rm -v "${map_dir}/${sample}.sam" \
#      "${map_dir}/${sample}.bam" \
      "${map_dir}/${sample}*nsort*"

EOF
done

# Hang script while there are still mapping jobs running
bmonitor "WGSBWA" 0
