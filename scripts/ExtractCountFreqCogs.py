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

    args = parser.parse_args()

#    import ipdb; ipdb.set_trace()

    contig_cog_list = defaultdict(list)
    for line in open(args.cog_file):
        line = line.rstrip()
        (contig, start, end, cog) = line.split(",")
        length = int(end) - int(start) + 1
        contig_cog_list[contig].append([int(start),int(end),int(length),cog])  
    
    
    for contig, cogs in contig_cog_list.iteritems():
        
        sorted_cogs = sorted(cogs, key=lambda tup: tup[0])
        
        contig_cog_list[contig] = sorted_cogs
    
    S = 0
    sample_files = []
    sample_names = []
    for name in glob.glob(args.input_dir + "/*.cnt.gz"):
        
        m = re.search(args.input_dir + "/(.*)\.cnt\.gz",name);
        if m:
            sample_name = m.groups()[0] 
            S = S + 1
            sample_files.append(name)
            sample_names.append(sample_name)
    
    for contig, sorted_cogs in contig_cog_list.iteritems():
        
        for cog in sorted_cogs:
            length = cog[2]
            
            freq_Array = np.zeros((length,S,4),dtype=np.int)
            
            cog.append(freq_Array)    


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
        
                    sorted_cogs = contig_cog_list[contig]
                
                    found = False
                    foundCog = None
                    for index, cog in enumerate(sorted_cogs):
                        if (pos >= cog[0] and pos <= cog[1]):
                            found = True
                            foundCog = cog
                            break 
                
                    if found:
                        for a in range(4): 
                            tok = ttokens[5 + a]

                            tokens2 = tok.split(":")
                            len2 = len(tokens2)
                            if len2 >= 14:
                                lpos = pos - foundCog[0] 
                                if lpos >= 0 and lpos < foundCog[2]:
                                    foundCog[4][lpos,s,a] = int(tokens2[1])
                                else:
                                    print "Cog index error: " + contig + " " + str(lpos)
                            else:
                                print "Format error: " + tok 
        s = s + 1
    

    expanded_sample_names = expand_sample_names(sample_names)

    results = []

    for contig, sorted_cogs in contig_cog_list.iteritems():
        
        for cog in sorted_cogs:
        
            length = cog[2]
        
            name = cog[3]
        
            cog_res = np.reshape(cog[4],(length,S*4))
        
            cog_df = p.DataFrame(cog_res,index=[name] * length)
        
            cog_df.columns = expanded_sample_names
            start = cog[0]
            #posns = [x + start for x in range(length)]
            posns = list(range(length))
            cog_df['Position'] = posns
            cols = cog_df.columns.tolist()
            cols = cols[-1:] + cols[:-1]
            cog_df = cog_df[cols]
        
            results.append(cog_df)
    
    result_df = p.concat(results)
    result_df.to_csv(args.output_file)

if __name__ == "__main__":
    main(sys.argv[1:])
