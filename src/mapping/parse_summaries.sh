# This script parses the custom summaries from bwa output into a structured
# table which can be used to generate plots with map_param.py
# Cyril Matthey-Doret
# 02.04.2017

echo 'alg,MM,K,W,tot,single,multi,miss' > mapstats.csv

for i in data/mapped/*;
do
    echo -n ${i##*/} | sed 's/\([a-zA-Z]*\)-\([0-9]*\)-\([0-9]*\)-\([0-9]*\)/\1,\2,\3,\4,/g' >> mapstats.csv
    awk 'BEGIN {FS = "\n"; OFS = ","; RS = ""; tot = 0; uniq = 0; multi = 0; miss = 0}
    { if ( NF==4 ) {
        { match($1,/: [0-9]+/); tot += substr($1, RSTART+2, RLENGTH) }
        { match($2,/: [0-9]+/); uniq += substr($2, RSTART+2, RLENGTH) }
        { match($3,/: [0-9]+/); multi += substr($3, RSTART+2, RLENGTH) }
        { match($4,/: [0-9]+/); miss += substr($4, RSTART+2, RLENGTH) }
    } }
    END {print tot, uniq, multi, miss}' $i/split_summary.log >> mapstats.csv
    
done
