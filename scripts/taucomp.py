

import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse
import cPickle

from operator import mul, div, eq, ne, add, ge, le, itemgetter
from itertools import izip
from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar

from scipy.stats import chi2

from sklearn.metrics import roc_curve, auc, accuracy_score


def variableTau(tau):
    """Calculates positions with variable bases"""
    N = tau.shape[0]
    G = tau.shape[1]
    variable_tau = np.zeros((N), dtype=bool)
    for v in range(N):
        diff = False
        id0 = np.where(tau[v,0,:] == 1)[0]
        for g in range(1,G):
            idg = np.where(tau[v,g,:] == 1)[0]
            if(idg[0] != id0[0]):
                diff = True 
            
            variable_tau[v] = diff
            
    return variable_tau

def compSND(tau1,tau2):
    G1 = tau1.shape[1]
    G2 = tau2.shape[1]
        
    snd = np.zeros((G1,G2),dtype=np.int)
    N = tau1.shape[0]
    for g in range(G1):
        #snd[g,g] = 0
        for h in range(G2):
            overlap = 0.0;
            for v in range(N):
                idg = np.where(tau1[v,g,:] == 1)[0]
                idh = np.where(tau2[v,h,:] == 1)[0]
                if(idg[0] == idh[0]):
                    overlap += 1 
                
            snd[g,h] = N - overlap
                
    return snd

def main(argv):

    #import ipdb; ipdb.set_trace()

    parser = argparse.ArgumentParser()
    
    parser.add_argument("gamma_file", help="relative frequencies of haplotypes")
    
    parser.add_argument("tau_file", help="haplotypes to determine uncertainty of")
        
    parser.add_argument('comp_files', type=file, nargs='+', 
                        help="tau_file2 ... tau_fileN replicates to compare to")
        
    args = parser.parse_args()
    
    gamma_file = args.gamma_file
    tau_file = args.tau_file
    comp_files = args.comp_files

    gamma = p.read_csv(gamma_file, header=0, index_col=0)
    gamma_mean = np.mean(gamma.as_matrix(),axis = 0) 

    tau = p.read_csv(tau_file, header=0, index_col=0)
    
    tau_matrix = tau.as_matrix()
    tau_matrix = np.delete(tau_matrix,0,1)
    
    V = tau_matrix.shape[0]
    G = tau_matrix.shape[1]/4
    
    tau_array = np.reshape(tau_matrix,(V, G,4)) 
    
    NC = len(comp_files)
    all_acc = np.zeros((G,NC))
    
    c = 0
    for ctau_file in comp_files:
    
        ctau = p.read_csv(ctau_file, header=0, index_col=0)
        ctau_matrix = ctau.as_matrix()
        ctau_matrix = np.delete(ctau_matrix,0,1)
    
        V2 = ctau_matrix.shape[0]
        G2 = ctau_matrix.shape[1]/4
    
        if V2 != V or G2 != G:
            print 'Haplotype files do not match V %d -> %d or G %d -> %d' % (V,V2,G,G2)
            sys.exit(-1)
    
        ctau_array = np.reshape(ctau_matrix,(V2, G2,4))
    
        comp = compSND(tau_array,ctau_array)/float(V)
    
        accuracies = np.zeros(G)
        map = np.zeros(G,dtype=int)
        acctotal = 0.0
        ga = 0;
        while (ga < G):
            (mr,mcol) = np.unravel_index(np.argmin(comp),comp.shape)
            curr_acc = np.min(comp)
            acctotal += curr_acc 
            comp[mr,:] = np.ones(G)
            comp[:,mcol] = np.ones(G)
            accuracies[mr] = curr_acc
            map[mr] = mcol
            ga += 1
        
        all_acc[:,c] = accuracies
        c=c+1
    
    mean_acc = np.mean(all_acc, axis=1)

    for g in range(G):
        print "%d,%f,%f\n" %(g,gamma_mean[g],mean_acc[g])
    
if __name__ == "__main__":
    main(sys.argv[1:])