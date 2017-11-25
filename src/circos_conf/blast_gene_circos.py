
# Extracting all BLAST hits between genes (non-self) and parse the results into
# A circos-compatible format.
# Cyril Matthey-Doret
# 24.11.2017
import pandas as pd
# Tabular output (fmt6) from BLAST
blast = "data/homology/MCScanX/input/MCScanX_in.blast"
# BED file with all genes coordinates
bed = "data/homology/MCScanX/input/MCScanX_genes_conv.bed"
# Significance threshold for E-values
Esig=10**-5
# Output folder
circos_links = 'data/circos/blast.lf.txt'
i=0
N = sum(1 for line in open(blast,'r'))

blast_file = pd.read_csv(blast, sep="\t", header=None)
bed_file = pd.read_csv(bed, sep="\t", header=None)
bed_file = bed_file.iloc[:,0:4]
signif = blast_file.loc[(blast_file[10] < Esig) & (blast_file[0] != blast_file[1])]
out_signif = signif.iloc[:,[0,1]]
Q_genes = bed_file.merge(out_signif, left_on = 3, right_on = 0)
QS_genes = bed_file.merge(Q_genes, left_on = 3, right_on = "1_y")
out_links = QS_genes.iloc[:,[1,2,3,5,6,7]]

out_links.to_csv(circos_links, sep = " ", header = False, index = False)
