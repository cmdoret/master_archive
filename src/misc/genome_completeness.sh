# Using busco to assess completeness of the genome
# Cyril Matthey-Doret
# 17.02.2018

source src/misc/dependencies.sh

# Prepare output folder
rm -rf data/ref_genome/busco/

# Retrieve database
wget -O data/ref_genome/busco_hymenoptera.tarhttp://busco.ezlab.org/datasets/hymenoptera_odb9.tar.gz
# run busco
run_BUSCO.py -i $REF -o busco -l hymenoptera_odb9.tar.gz -m geno -c 5

# Clean up
rm hymenoptera_odb9.tar.gz
mv busco* data/ref_genome/busco/
