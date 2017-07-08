#!bin/bash

# This script parses the output of pstacks into a table to determine the best value for minimum coverage.
# Extracted metrics are: Mean coverage, st dev coverage, number of loci
# Cyril Matthey-Doret
# 11.04.2017



cd "$(dirname "$0")"

pdir=../../data/logs/pstacks
echo "min_cov,nloci,mean_cov,sd_cov" > pstats.csv

for f in $pdir/PST*;
do
    echo -n $(echo $f | sed 's/.*\([0-9]\).*/\1/')"," >> pstats.csv
    awk 'BEGIN {FS = "\n"; OFS = ","; RS = ""; rds = 0; meancov = 0; sdcov = 0; nloc = 0; nstacks = 0; n = 0}
    {
    { match($5,/Analyzed [0-9]+/);AR = substr($5,RSTART+9,RLENGTH-9) }
    if ( AR > 0 ) {
        { match($11,/Wrote [0-9]+/);nloc += substr($6,RSTART+6,RLENGTH-6) }
        { match($9,/is [0-9]+/);meancov += substr($9,RSTART+3,RLENGTH-3) }
        { match($9,/Dev: [0-9]+/);sdcov += substr($9,RSTART+5,RLENGTH-5) }
        {n+=1}
    }}
    END {printf "%.1f,%.1f,%.1f\n", nloc/n,meancov/n,sdcov/n}' $f >> pstats.csv
done;
cp pstats.csv ../../reports/lab_book/
