# Generate QC reports for both fastq files and Mapping.
# Needs to be run after the mapping step.
# Intended for the analysis of wild samples WGS data.
# Cyril Matthey-Doret
# 23.02.2018

# Parse CL arguments
while [[ "$#" > 1 ]];
do
    case $1 in
        # Working directory. Must contain a "raw" folder with reads
        --workdir) wgs="$2";;
        # Reference genome. Must be indexed beforehand
        --ref) index="$2";;
        # Output folder
        --out) out_folder="$2";;
        *) break;;
esac; shift; shift
done

source src/misc/jobs_manager.sh
MAXPROC=15

# Folder containing input raw reads files
in_dir="${wgs}/raw/"
map_dir="${wgs}/mapped/"

# Resetting temporary working directory
tmp_dir="${wgs}/tmp/"
rm -rf "$tmp_dir"
mkdir -p "$tmp_dir"

logs="${wgs}/log/"
mkdir -p "$logs"


# List of sample names
samples=( $(find "$in_dir" | grep "fastq" | sed 's/.*\(HYI-[0-9]*\)_R.*/\1/g' | sort | uniq) )

for sample in ${samples[@]};
do

bmonitor "WGSQC" $MAXPROC

# Submit the QC of each sample as an independent job
bsub << EOF
#!/bin/bash
#BSUB -L /bin/bash
#BSUB -o ${logs}/QC-${sample}-OUT.txt
#BSUB -e ${logs}/QC-${sample}-ERROR.txt
#BSUB -u cmatthey@unil.ch
#BSUB -J WGSQC-${sample}
#BSUB -n 4
#BSUB -R "span[ptile=4]"
#BSUB -q priority
#BSUB -R "rusage[mem=8000]"
#BSUB -M 8000000

# Loading softwares
source src/misc/dependencies.sh

sample_fq="${tmp_dir}/${sample}.fastq.gz"

# Merge all fastq files for each sample (both lanes, reverse and forward)
find "${in_dir}" -name "*${sample}*" -type f | sort | xargs cat > "\$sample_fq"

# Fastqc on each sample
fastqc \$sample_fq -o $out_folder

# mapping stats with samtools stats
samtools stats -r $index "${map_dir}/${sample}.fixed.csort.bam" > $out_folder/${sample}_stat.txt

EOF
done

bmonitor "WGSQC" 0

# MultiQC for all samples: Incorporating samtools picard and fastqc output
multiqc ${out_folder}/* ${map_dir}/*dup_metrics* -o $out_folder

# remove temporary fq
rm -rf "$tmp_dir"
