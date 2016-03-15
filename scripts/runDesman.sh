#!/bin/bash

varFile=$1

contig_6or16_genesL_sel_var.csv -e contig_6or16_genesL_tran_df.csv


eFile=$2
oFile=$3

somecommand ${1:-foo}

for g in 2 3 4 5 6 7 8  
do
	for r in 0 1 2 3 4
        do
        	desman $varFile -e $eFile -o contig_6or16_${g}_${r} -r 1000 -i 100 -g $g -s $r > contig_6or16_${g}_${r}.out&  
        done
done
