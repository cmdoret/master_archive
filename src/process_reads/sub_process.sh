#!bin/bash

LIB=$1
case $# in 
1)
MM=2
;;
2)
MM=$2
;;
esac


bsub <<EOF
$(bash 'process_reads/process_'$LIB'.sh' $MM)
EOF
