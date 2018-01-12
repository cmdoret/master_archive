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
        --workdir) WGS="$2";;
        # Reference genome. Must be indexed beforehand
        --ref) IDX="$2";;
        *) break;;
esac; shift; shift
done

# #Mismatches allowed per read for mapping
MM=4
# #threads used when mapping
threads=4
# Max number of parallel jobs
MAXPROC=5

# Allows to cap # parallel jobs and wait for their execution
source src/misc/jobs_manager.sh


wgs="$WGS"
index="$IDX"

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

declare -i progress=0
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
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -q normal
#BSUB -R "rusage[mem=64000]"
#BSUB -M 64000000


module add UHTS/Aligner/bwa/0.7.13;
module add UHTS/Analysis/samtools/1.3;

forward="${wgs}/merged/${sample}_R1.fastq.gz"
reverse="${wgs}/merged/${sample}_R2.fastq.gz"

# Merge lanes, keep forward and reverse separated
find "${in_dir}" -name "*${sample}*R1*" -type f | xargs cat > "\$forward"
find "${in_dir}" -name "*${sample}*R2*" -type f | xargs cat > "\$reverse"

# Mapping paired ends separately
bwa aln -n $MM -t $threads "$index" \$forward > "${tmp_dir}/${sample}_R1.sai"
bwa aln -n $MM -t $threads "$index" \$reverse > "${tmp_dir}/${sample}_R2.sai"

# Combine information from forward and reverse to generate sam files
bwa sampe $index "${tmp_dir}/${sample}_R1.sai" "${tmp_dir}/${sample}_R2.sai" \
          "\$forward" "\$reverse" > "${map_dir}/${sample}.sam"

# Splitting reads into not mapped, unique map and multi map
perl src/mapping/split_sam.pl -i "${map_dir}/${sample}.sam" -o "${map_dir}/${sample}" >> "${map_dir}/split_summary.log"

# Removing original sam file
rm -v "${map_dir}/${sample}.sam"

# Convert SAM files to BAM
mkdir -p "${map_dir}/bam/"
samtools view -@ $threads -bS -o "${map_dir}/bam/${sample}-uniq.bam" "${map_dir}/${sample}-uniq.sam"

# Sort alignments by leftmost coordinate
samtools sort -@ $threads "${map_dir}/bam/${sample}-uniq.bam" -o "${map_dir}/bam/${sample}-uniq.sorted.bam"

# Index BAM files
samtools index "${map_dir}/bam/${sample}-uniq.sorted.bam"

# Output index statistics
samtools idxstats "${map_dir}/bam/${sample}-uniq.sorted.bam"

# Remove unsorted bam files
rm -v "${map_dir}/bam/${sample}-uniq.bam"
# Compress all SAM files
gzip -v "${map_dir}/${sample}*.sam"
EOF

((progress++))
done

# Hang script while there are still mapping jobs running
bmonitor "WGSBWA" 0
