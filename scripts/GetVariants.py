import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss

from scipy.stats import norm

import argparse
import math
import re

from Bio import SeqIO
from collections import defaultdict

def read_gene_file(gene_file): 
    
    gene_mapping = {}
    #contig-1001153_3,COG4577,2749,3034,2.0,90.0,1
    for line in open(gene_file):
        line = line.rstrip()
    
        (gene_id,cog,start,end,a,b,strand) = line.split(",")
    
        m = re.search(r"(.*)_\d+", gene_id)        
        contig = m.group(1)
        gene_mapping[cog]= (contig,int(start),int(end),int(strand))
        
    return gene_mapping

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("contig_fasta_file", help="contig sequences")

    parser.add_argument("gene_file", help="gene locations")
    
    parser.add_argument("tau_file", help="variants")

    args = parser.parse_args()

    import ipdb; ipdb.set_trace()

    base_map = ['A','C','G','T']
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

    gene_mapping = read_gene_file(args.gene_file)

    contigs = {}

    handle = open(args.contig_fasta_file, "rU")

    for record in SeqIO.parse(handle, "fasta"):
	    seq = record.seq
	    contigs[record.id] = str(seq)
    
    tau_df = p.read_csv(args.tau_file, header=0, index_col=0)
    
    pos = tau_df['Position']
    gene_ids = set(tau_df.index.tolist())
    
    G = (len(tau_df.columns) - 1)/4
    c = 0
    C = len(gene_ids)
    
    variants = defaultdict(dict)
    for gene in gene_ids:
        
        try:
            gene_variants = tau_df.loc[gene]
                
            if len(gene_variants.shape) == 1:
                gene_variants = gene_variants.to_frame()
                gene_variants = gene_variants.transpose()
        
            gene_pos = gene_variants['Position']
            gene_variants_matrix = gene_variants.as_matrix()
            gene_variants_matrix = np.delete(gene_variants_matrix, 0,1)
            V = gene_variants_matrix.shape[0]
            gene_variants_matrix = np.reshape(gene_variants_matrix,(V,G,4))
        
            for g in range(G):
                (contig,start,end,strand) = gene_mapping[gene]
                contig_seq = contigs[contig]
                variants[gene][g] = list(contig_seq[start:end])
            v = 0
            for (index,value) in gene_pos.iteritems():
                value = value - 1
                for var in variants[gene]:
                    base = np.nonzero(gene_variants_matrix[v,int(var)])[0]
                  
                    variants[gene][var][value] = base_map[base[0]]
                v = v + 1
        except KeyError:
            for g in range(G):
                (contig,start,end,strand) = gene_mapping[gene]
                contig_seq = contigs[contig]
                variants[gene][g] = list(contig_seq[start:end])
            pass
        
        c = c + 1
    
    for gene in gene_ids:
        (contig,start,end,strand) = gene_mapping[gene]
        
        name = "Cogs/" + gene + "_hap.fa"
        
        f = open(name,'w')
        
        for var in variants[gene]:
            seq_list = variants[gene][var]
            
            if strand > 0:
                seq = ''.join(seq_list)
            else:
                comp_list = [complement[base] for base in seq_list] 
                comp_seq = ''.join(comp_list)
                seq = comp_seq[::-1]
            
            f.write(">" + "Haplotype" + "_" + str(var) + "\n")
            f.write(seq + "\n")
            
if __name__ == "__main__":
    main(sys.argv[1:])
