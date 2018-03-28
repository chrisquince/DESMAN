import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
import argparse
import math
import logging

from operator import mul, truediv as div, eq, ne, add, ge, le, itemgetter
from numpy.random import RandomState


def compGenes(etaPred, etaG):
    G = etaPred.shape[1]
    NC = etaPred.shape[0]
    NG = etaG.shape[1]
    
    acc = np.zeros((NG,G))
    for g in range(NG):
        for h in range(G):
            acc[g,h] = np.sum(etaG[:,g] == etaPred[:,h])
            
            acc[g,h] =acc[g,h]/float(NC)
    accR = np.copy(acc)            
    accuracies = np.zeros(NG)
    map = np.zeros(NG,dtype=int)
    acctotal = 0.0
    ga = 0;
    while (ga < NG):
        (mr,mcol) = np.unravel_index(np.argmax(acc),acc.shape)
        curr_acc = np.max(acc)
        acctotal += curr_acc 
        acc[mr,:] = np.zeros(G)
        acc[:,mcol] = np.zeros(NG)
        accuracies[mr] = curr_acc
        map[mr] = mcol
        ga += 1
    
    acctotal = acctotal/float(NG)
    return (acctotal,accuracies,accR)


def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("assign_file", help="input core gene coverages")

    parser.add_argument("genomes", help="input MAP estimate frequencies")

    args = parser.parse_args()   
   
    #import ipdb; ipdb.set_trace()

    genomes    = p.read_csv(args.genomes, header=0, index_col=0)
    assigns    = p.read_csv(args.assign_file, header=0, index_col=0)
    eta     = assigns.as_matrix()
    contig_names = assigns.index.tolist()

    genomes = genomes.loc[contig_names]
    genomes_M   = genomes.as_matrix()
    genomes_D = np.copy(genomes_M)
        
    genomes_D[genomes_D < 0.5] = 0.
    genomes_D[genomes_D >= 0.5] = 1.0
        
    eta[eta < 0.5] = 0.
    eta[eta >= 0.5] = 1.0
        
    (dtotal, dacc,dacc_array) = compGenes(eta, genomes_D)
    
    np.savetxt(sys.stdout, dacc, fmt='%.4f')
    
    print(('Av. accurracy = %f' %(dtotal)))
       
if __name__ == "__main__":
    main(sys.argv[1:])
