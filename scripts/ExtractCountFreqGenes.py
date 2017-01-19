import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
import argparse
import math
import glob
import re
import gzip

from Bio import SeqIO

from collections import defaultdict
from collections import Counter
from scipy import stats

def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring

def expand_sample_names(sample_names):
    expanded = []
    
    for name in sample_names:
        expanded.append(name + "-A")
        expanded.append(name + "-C")
        expanded.append(name + "-G")
        expanded.append(name + "-T")
    
    return expanded

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("cog_file", help="cogs and contig locations for frequencies to be called on")

    parser.add_argument("input_dir", help="input directory to glob *.cnt.gz from")

    parser.add_argument("--output_file", default="Select_freq.csv")
    
    parser.add_argument('-g', '--gene_file', action='store_true',help=("alternate input format"))
    
    args = parser.parse_args()
    

    #import ipdb; ipdb.set_trace()

    contig_gene_list = defaultdict(list)
    for line in open(args.cog_file):
        line = line.rstrip()
    #COG0100,k99_17880,11324,11714,k99_17880_25,1

    #k99_10285_3,k99_10285,873,1386,1
        if args.gene_file:
            (gene, contig, start, end, strand) = line.split(",")
        else:
            (cog, contig, start, end, gene, strand) = line.split(",")

        length = int(end) - int(start) + 1
        contig_gene_list[contig].append([int(start),int(end),int(length),gene])  
    
    
    for contig, genes in contig_gene_list.iteritems():
        
        sorted_genes = sorted(genes, key=lambda tup: tup[0])
        
        contig_gene_list[contig] = sorted_genes
    
    S = 0
    sample_files = []
    sample_names = []
    abs_path_input_dir = os.path.abspath(args.input_dir)
    
    for name in glob.glob(abs_path_input_dir + "/*.cnt.gz"):
        basename = os.path.basename(name)    
        m = re.search("(.*)\.cnt\.gz",basename)
        if m:
            sample_name = m.groups()[0] 
            S = S + 1
            sample_files.append(name)
            sample_names.append(sample_name)
    
    for contig, sorted_genes in contig_gene_list.iteritems():
        
        for gene in sorted_genes:
            length = gene[2]
            
            freq_Array = np.zeros((length,S,4),dtype=np.int)
            
            gene.append(freq_Array)    


    s = 0
    for file_name in sample_files:
        
        print "Reading frequencies from " + file_name
        
        with gzip.open(file_name,'r') as fin:
            for line in fin:
                line = line.rstrip()
                
                ttokens = line.split('\t')
                ltokens = len(ttokens)
                if ltokens < 10:
                    print "Format error: " + line
                else:
                    contig = ttokens[0]
                    pos = int(ttokens[1])
        
                    sorted_genes = contig_gene_list[contig]
                
                    found = False
                    foundGene = None
                    for index, gene in enumerate(sorted_genes):
                        if (pos >= gene[0] and pos <= gene[1]):
                            found = True
                            foundGene = gene
                            break 
                
                    if found:
                        for a in range(4): 
                            tok = ttokens[5 + a]

                            tokens2 = tok.split(":")
                            len2 = len(tokens2)
                            if len2 >= 14:
                                lpos = pos - foundGene[0] 
                                if lpos >= 0 and lpos < foundGene[2]:
                                    foundGene[4][lpos,s,a] = int(tokens2[1])
                                else:
                                    print "Cog index error: " + contig + " " + str(lpos)
                            else:
                                print "Format error: " + tok 
        s = s + 1
    

    expanded_sample_names = expand_sample_names(sample_names)

    results = []

    for contig, sorted_genes in contig_gene_list.iteritems():
        
        for gene in sorted_genes:
        
            length = gene[2]
        
            name = gene[3]
        
            gene_res = np.reshape(gene[4],(length,S*4))
        
            gene_df = p.DataFrame(gene_res,index=[name] * length)
        
            gene_df.columns = expanded_sample_names
            start = gene[0]
            posns = list(range(length))
            gene_df['Position'] = posns
            cols = gene_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            gene_df = gene_df[cols]
        
            results.append(gene_df)
    
    result_df = p.concat(results)
    result_df.to_csv(args.output_file)

if __name__ == "__main__":
    main(sys.argv[1:])
