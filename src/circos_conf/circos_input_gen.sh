# This script generates the input files required to run circos.
# Cyril Matthey-Doret
# 21.11.2017
#---SET UP DIRECTORY---#
ref="$1"
gwas="data/assoc_mapping/case_control/case_control_all.tsv"
cir_dir="$2"
rm -rf "$cir_dir"
mkdir -p "$cir_dir"

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

# GWAS Manhattan plot
# parsing p-values of SNPS into: chr start end value
# where start=end (points)
tail -n +2 $gwas |
  awk 'BEGIN{IFS="\t";OFS=" "}
       {print $2,$3,$3,$13}' |
  grep "chr" |
  sed 's/chr/lf/' > "$cir_dir/gwas.lf.txt"

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
<plots>
<<include $cir_dir/lf.centro.conf>>
<<include $cir_dir/lf.gwas.conf>>
</plots>
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
thickness        = 20p
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

glyph            = rectangle
glyph_size       = 8
color            = grey
stroke_color     = grey
stroke_thickness = 1

<axes>
<axis>
color = lred
thickness = 3
position = 3
</axis>
</axes>


<backgrounds>
<background>
color     = greys-3-seq-1
y0 = 0
</background>
</backgrounds>

<rules>
<rule>
condition    = var(value) > 3
color        = dred
fill_color   = dred_a1
glyph_size       = 14
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

#<link>
#file          = data/5/segdup.txt
#radius        = 0.8r
#bezier_radius = 0r
#color         = black_a4
#thickness     = 2
#</link>

</links>
MCSCANX

# OPTIONAL
# contig boundaries
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
