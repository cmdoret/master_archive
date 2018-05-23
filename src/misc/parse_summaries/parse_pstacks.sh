#!bin/bash

# This script parses the output of pstacks into a table to determine the best value for minimum coverage.
# Extracted metrics are: Mean coverage, st dev coverage, number of loci
# Cyril Matthey-Doret
# 11.04.2017


# Working directory set to script location
cd "$(dirname "$0")"

# Folder containing pstacks log files
pdir=../../data/logs/pstacks

# Iterating over range of coverage parameter values
for f in {1..6};
do
    # If the lofile exist for this parameter value
    if [ $(ls -l $pdir | grep "COVMIN${f}" | wc -c) -gt 0 ];
    then
        # Sending content from all samples' STDERR to single file
        awk 'FNR==1{print ""}{print}' $pdir/PST_COVMIN${f}*ERR* > $pdir/summary_PST${f}.txt  || echo "no pstacks log files for m=$f"
    fi;
done;

# Writing header of summary table output file
echo "min_cov,nloci,mean_cov,sd_cov" > pstats.csv

for f in $pdir/summary*;
do
    # For each parameter value tested initiate line with parameter value 
    echo -n $(echo $f | sed 's/.*\([0-9]\).*/\1/')"," >> pstats.csv
    Send the
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
cp pstats.csv ../../docs/lab_book/
