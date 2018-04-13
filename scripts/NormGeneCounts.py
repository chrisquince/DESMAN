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

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gene_counts_file", help="input gene base frequencies")

    parser.add_argument('-f','--frac_thresh', type=float, default=0.01,
                 help=("threshold for assigning gene to strain as fraction of reads deriving from it, defaults to 0.01"))

    parser.add_argument("--output_file", default="Gene_etaP.csv")

    args = parser.parse_args()

    import ipdb; ipdb.set_trace()

    gene_counts    = p.read_csv(args.gene_counts_file, header=0, index_col=0, sep='\t')
    gene_counts_matrix = (gene_counts.as_matrix()).astype(float)
    row_sums = gene_counts.sum(axis=1)
    gene_countsP = gene_counts_matrix / row_sums[:, np.newaxis]
    
    gene_countsP[gene_countsP > args.frac_thresh] = 1.
    gene_countsP[gene_countsP <= args.frac_thresh] = 0.

    gene_countsP_df = p.DataFrame(gene_countsP,index=gene_counts.index.get_values())

    gene_countsP_df.to_csv(args.output_file)

if __name__ == "__main__":
    main(sys.argv[1:])
