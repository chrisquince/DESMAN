#!/bin/bash

varFile=$1

stub=${varFile%sel_var.csv}

eFile=${2:-${stub}tran_df.csv}
oFile=${3:-$stub}

echo $stub
echo $eFile
echo $oFile

for g in 2 3 4 5 6 7 8  
do
	for r in 0 1 2 3 4
        do
		echo $g
        	desman $varFile -e $eFile -o ${stub}_${g}_${r} -i 500 -g $g -s $r -r 1000 > ${stub}_${g}_${r}.out&  
        done
done
