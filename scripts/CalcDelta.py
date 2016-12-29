#!/usr/bin/env python
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
    
    originalS = (len(scgCols) - 1)/4
    
    idx = range(1,originalS*4,4)
    
    sampleNames = [scgCols[i] for i in idx] 

    return sampleNames

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gene_cov_file", help="input gene base frequencies")

    parser.add_argument("gene_scg_file", help="core genes")

    parser.add_argument("output_stub", help="output file stub")

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()

    scgs = []

    for line in open(args.gene_scg_file):
        line = line.rstrip()
        scgs.append(line)
    
    scgs = set(scgs)
    NS = len(scgs)
    scgs_covs = None 
    s = 0
    header = None
    for line in open(args.gene_cov_file):
        line = line.rstrip()
        
        if header is None:
            header = line
            samples = line.split(",")
            samples.pop(0)
            scgs_covs = np.zeros((NS,len(samples)))
        else:
            covs = line.split(",")
    
            gene = covs.pop(0)
            if gene in scgs:
                scgs_covs[s,:] = covs 
                s = s+1
    total_mean = np.mean(scgs_covs,axis=0,dtype=np.float64)
    total_sd = np.std(scgs_covs,axis=0,dtype=np.float64)
    
    mean_sd =  np.column_stack((total_mean,total_sd))
    mean_sd_df = p.DataFrame(mean_sd,index=samples)
    mean_sd_df.columns = ['mean','sd']
    mean_sd_df.to_csv(args.output_stub+"mean_sd_df.csv")
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(total_mean,total_sd)
    
if __name__ == "__main__":
    main(sys.argv[1:])
