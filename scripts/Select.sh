#!/bin/bash

while IFS=, read cog contig start end gene 
do
    cp Hits/${cog}.gfa Select
    
done < ../Annotate/ClusterEC_core.cogs 

