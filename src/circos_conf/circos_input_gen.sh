# This script generates the input files required to run circos.
# Cyril Matthey-Doret
# 21.11.2017

#---SET UP DIRECTORY---#
ref="data/ref_genome/ordered_genome/merged.fasta"
tig="data/annotations/corresp_gff.csv"
gwas="data/assoc_mapping/case_control/case_control_all.tsv"
gwas_hit="data/assoc_mapping/case_control/case_control_hits.tsv"
mcsx="data/homology/MCScanX/input/MCScanX_in.collinearity"
gff="data/homology/MCScanX/input/MCScanX_genes_conv.gff"
cir_dir="data/circos/"
rm -rf "$cir_dir"
mkdir -p "$cir_dir"
# power of 10 to consider hits significant
sigpow=2

#---FORMAT DATA FILES---#

# Karyotype file:
# Parsing FASTA reference to karyotype format: chr - name label start end color
awk '$0 ~ ">" {
         print c; c=0;printf substr($0,2,100) "\t"; }
     $0 !~ ">" {
         c+=length($0);}
       END { print c; }' $ref |
  grep "chr" |
  sed 's/chr/lf/' |
  awk 'BEGIN{
         OFS=" ";N=0}
       {
         N+=1;print "chr","-",$1,N,0,$2,"chr"N}' > "$cir_dir/karyotype.lf.txt"

# Contig boundaries
tail -n +2 $tig |
  awk 'BEGIN{
         FS=",";OFS=" "}
       $2 ~ "chr" {
         if ( $5 ~ "rev" ){print $2,$3-$4,$3}
       else {print $2,$3,$3+$4}
       }' |
  sed 's/chr/lf/' |
  sort -k1,1 -k2,2n > "$cir_dir/tig_lim.lf.txt"

# Contigs with CSD hits
while read line
do
  # is the p-value significant ?
  pval=$(echo $line | awk '{print $NF}' | sed 's/[eE]+\{0,1\}/*10^/g')
  if (( $(echo "$pval  > $sigpow" | bc -l) ))
  then
    awk -v chrom="$(cut -d$' ' -f2 <(echo $line))" \
        -v pos="$(cut -d$' ' -f3 <(echo $line))" \
        '$1 == chrom {
          if ( $2 <= pos && $3 >= pos ) {print $0}}' "$cir_dir/tig_lim.lf.txt"
  fi
done < <( tail -n +2 "$gwas_hit" | sed 's/chr/lf/') |
  uniq > "$cir_dir/csd_tig.lf.txt"

# GWAS Manhattan plot
# parsing p-values of SNPS into: chr start end value
# where start=end (points)
tail -n +2 $gwas |
  awk 'BEGIN{IFS="\t";OFS=" "}
       {print $2,$3,$3,$13}' |
  grep "chr" |
  sed 's/chr/lf/' > "$cir_dir/gwas.lf.txt"

# Collinearity blocks
# Parsing MCScanX collinearity output into circos links
while read line
do
  awk -v g1="${line%% *}" -v g2="${line##* }" \
    'BEGIN{FS=" ";N1=0;N2=0}
     $2 ~ g1 {
       N1+=1;c1=$1;s1=$3;e1=$4}
     $2 ~ g2 {
       N2+=1;c2=$1;s2=$3;e2=$4}
     END{
       if (N1 != 1 || N2 != 1) {exit 1}
       else {print c1,s1,e1,c2,s2,e2}}' $gff
  # Convert gene IDs to coordinates
done < <(grep -v "^#" $mcsx | tr -d ' ' | tr '\t' ' ' | cut -d$' ' -f2,3) |
  sed 's/chr/lf/g' > "$cir_dir/mcsx.lf.txt"

# Filtering links in/outside CSD contigs
echo -n "" > "$cir_dir/mcsx_csd.lf.txt"
echo -n "" > "$cir_dir/mcsx_out.lf.txt"
while read -a line
do
  echo $lineawk
done < "$cir_dir/mcsx.lf.txt"


#---GENERATE CONFIGURATION---#
# config file:
cat << CFG > $cir_dir/lf.main.conf
<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<<include $cir_dir/lf.ideogram.conf>>
<<include $cir_dir/lf.ticks.conf>>
<<include $cir_dir/lf.mcscanx.conf>>
karyotype = "$cir_dir/karyotype.lf.txt"
chromosomes_units = 1000000
# Standard stuff:
<image>
<<include etc/image.conf>>
</image>
CFG

cat << IDEOGRAM > $cir_dir/lf.ideogram.conf
<ideogram>
<spacing>
# Spacing between ideograms. Suffix "r" denotes a relative value. It
# is relative to circle circumference (e.g. space is 0.5% of
# circumference).
default = 0.005r
#<pairwise lf1;lf2>
#spacing = 20r
#</pairwise>

</spacing>
radius           = 0.90r
thickness        = 60p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p

# LABELS
show_label       = yes
label_font       = default
label_radius     = 1r + 75p
label_size       = 50
label_parallel   = yes

</ideogram>

<highlights>

z          = 5

 <highlight>
 # Contig boundaries
 ideogram   = yes
 file       = "$cir_dir/tig_lim.lf.txt"
 fill_color = blue
 stroke_color = black
stroke_thickness = 2
 </highlight>

 <highlight>
 # CSD contig
 ideogram   = yes
 file       = "$cir_dir/csd_tig.lf.txt"
 z          = 5
 fill_color = red
 </highlight>
</highlights>

<plots>
<plot>
# CSD highlight
type = highlight
fill_color = red_a5
stroke_color = red_a5
file       = "$cir_dir/csd_tig.lf.txt"
r0   = 0.80r
r1   = 1r
z    = 10
</plot>
<<include $cir_dir/lf.centro.conf>>
<<include $cir_dir/lf.gwas.conf>>
</plots>

IDEOGRAM

cat << GWAS > $cir_dir/lf.gwas.conf
#GWAS p-values (scatter plot)
<plot>

show  = yes
type  = scatter

file  = $cir_dir/gwas.lf.txt
r1    = 0.98r
r0    = 0.80r
max   = 6.2.0
min   = 0.0
orientation = in

glyph            = circle
glyph_size       = 8
color            = lgrey
stroke_color     = lgrey
stroke_thickness = 1

<axes>
<axis>
color = grey
thickness = 3
position = $sigpow
</axis>
</axes>

<backgrounds>
<background>
color = prgn-3-div-2
y0 = 0
</background>
</backgrounds>

<rules>
<rule>
condition    = var(value) > $sigpow
color        = dred
fill_color   = dred_a1
glyph_size       = 18
</rule>
</rules>

</plot>

GWAS

cat << CENTRO > $cir_dir/lf.centro.conf
# Centromere loess (line)
CENTRO

cat << MCSCANX > $cir_dir/lf.mcscanx.conf
# Collinearity blocks
<links>

<link>
file          = "$cir_dir/mcsx.lf.txt"
radius        = 0.8r
bezier_radius = 0r
color         = dgrey_a4
thickness     = 2
</link>

</links>
MCSCANX

# OPTIONAL
# coverage / FPKM ?
# allelic diversity
# transcripts



cat << TICKS > $cir_dir/lf.ticks.conf
# TICKS

show_ticks          = yes
show_tick_labels    = yes

<ticks>

radius           = 1r
color            = black
thickness        = 2p

multiplier       = 1e-6

format           = %d

<tick>
show_label = yes
label_size = 25
spacing        = 5u
size           = 20p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>
</ticks>
TICKS

#---RUN CIRCOS---#
export PATH="~/Public/circos-0.69-6/bin/":$PATH
circos -conf "$cir_dir/lf.main.conf"
