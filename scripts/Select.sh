#!/bin/bash

while IFS=, read cog contig start end gene 
do
    cp Hits/${cog}.gfa Select
done < ../AnnotateEC/ClusterEC_core.cogs 

