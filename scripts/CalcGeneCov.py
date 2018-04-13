#!/usr/bin/env python3


import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
import argparse
import math

from collections import defaultdict
from collections import Counter
from scipy import stats


def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring

def get_sample_names(scg_freq):
    scgCols = scg_freq.columns.values.tolist()
    
    originalS = (len(scgCols) - 1) // 4
    
    sampleNames = [scgCols[i] for i in range(1,originalS*4,4)] 

    return sampleNames

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gene_freq_file", help="input gene base frequencies")

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()

    #read in data
    gene_freq  = p.read_csv(args.gene_freq_file, header=0, index_col=0)
    
    gene_set = set(gene_freq.index.tolist())
    genes = list(gene_set)
    N_genes = len(genes)
    gene_covs = {}
    gene_N = {}
    gene_stds = {}
    c = 0
    
    sample_names = get_sample_names(gene_freq)
    sample_names_new = [ rchop(x,'-A') for x in sample_names ] 
    
    S = len(sample_names)
    
    all_covs = np.zeros(N_genes*S) 
    all_stds = np.zeros(N_genes*S)
    
    g = 0   
    for gene in genes:
        try:
            gene_freqs = gene_freq.loc[gene]
            gene_freq_matrix = gene_freqs.as_matrix()
            gene_freq_matrix = np.delete(gene_freq_matrix, 0, 1)
            gene_freqs = np.reshape(gene_freq_matrix, (gene_freq_matrix.shape[0],gene_freq_matrix.shape[1] // 4,4))
            gene_totals = gene_freqs.sum(axis=2)
    
            total_mean = np.mean(gene_totals,axis=0,dtype=np.float64)
            total_sd = np.std(gene_totals,axis=0,dtype=np.float64)
            
            gene_covs[gene] = total_mean
            gene_stds[gene] = total_sd
            
            
            all_covs[g*S:(g+1)*S] = total_mean
            all_stds[g*S:(g+1)*S] = total_sd
            g = g + 1   
            pass
        except KeyError:
            gene_covs[gene] = None
            gene_stds[gene] = None
            gene_N[gene] = 0

            pass
    
    #slope, intercept, r_value, p_value, std_err = stats.linregress(all_covs,all_stds)
    
    sString = ",".join(sample_names_new)
    print("Gene," + sString)
    for gene in genes:
        covs = gene_covs[gene]    
        cString = ','.join(map(str, covs)) 

        print(gene + "," + cString)
        
   
    
    
    
if __name__ == "__main__":
    main(sys.argv[1:])
