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
    
    originalS = (len(scgCols) - 1)/4
    
    idx = list(range(1,originalS*4,4))
    
    sampleNames = [scgCols[i] for i in idx] 

    return sampleNames

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gene_freq_file", help="input gene base frequencies")

    parser.add_argument("select_cov_file", help="input gene base frequencies")

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()

    #import ipdb; ipdb.set_trace()

    select = []

    header1 = None
    for line in open(args.select_cov_file):
        line = line.rstrip()
        
        if header1 is None:
            header1 = line
        else:
            samples = line.split(",")
            gene = samples.pop(0)
            select.append(gene)
    
    select = set(select)
    
    header = None
    for line in open(args.gene_freq_file):
        line = line.rstrip()
        
        if header is None:
            header = line
            print(line)
            
        else:
            covs = line.split(",")
    
            gene = covs.pop(0)
            if gene in select:
                print(line)
        
   
    
    
    
if __name__ == "__main__":
    main(sys.argv[1:])