#!/usr/bin/env python3
import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
import argparse
import math

from Bio import SeqIO

from collections import defaultdict
from collections import Counter
from scipy import stats

def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring



def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("coverage_file", help="contig coverage csv")

    parser.add_argument("cluster_file", help="cluster assignments")
    
    parser.add_argument("fasta_file", help="fasta sequences")

    args = parser.parse_args()

#    import ipdb; ipdb.set_trace()

    #read in coverages
    coverage  = p.read_csv(args.coverage_file, header=0, index_col=0)
    coverage_matrix = coverage.as_matrix()
    
    contig_names = coverage.index.values
    sample_names = list(coverage.columns.values)        
    contig_dict = dict(list(zip(contig_names,list(range(len(contig_names))))))
        
    cluster_mapping = {}
    
    for line in open(args.cluster_file):
        line = line.rstrip()
    
        tokens = line.split(",")
    
        cluster_mapping[tokens[0]] = tokens[1]
        
    
    contig_lengths = {}
    handle = open(args.fasta_file, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        seq = record.seq
        
        contig_lengths[record.id] = len(seq)
    handle.close()
    
    #sum coverages in each cluster
    NSamples = coverage_matrix.shape[1]
    clusterSum = defaultdict(lambda: np.zeros(NSamples))
    clusterLengthSum = Counter()
    for contig, cluster in cluster_mapping.items():
        idx = contig_dict[contig]
        clusterSum[cluster] += coverage_matrix[idx,:]*contig_lengths[contig]
        clusterLengthSum[cluster] += contig_lengths[contig]

    hString = ",".join(sample_names)
    print("Cluster" + "," + hString)
    for cluster, sum in clusterSum.items():
        sumR = sum/clusterLengthSum[cluster]
        sumString = ",".join(map(str,sumR.tolist()))
        print(cluster + "," + sumString)
    

    #import ipdb; ipdb.set_trace()
if __name__ == "__main__":
    main(sys.argv[1:])
