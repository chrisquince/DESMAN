

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
    parser.add_argument("tau_file1", help="predicted variants")
        
    parser.add_argument("tau_file2", help="known variants")
        
    args = parser.parse_args()
    
    tau_file1 = args.tau_file1
    tau_file2 = args.tau_file2

    tau1 = p.read_csv(tau_file1, header=0, index_col=0)
    
    tau1_matrix = tau1.as_matrix()
    tau1_matrix = np.delete(tau1_matrix,0,1)
    
    V1 = tau1_matrix.shape[0]
    G1 = tau1_matrix.shape[1]/4
    
    tau1_array = np.reshape(tau1_matrix,(V1, G1,4)) 
    
    tau2 = p.read_csv(tau_file2, header=0, index_col=0)
    tau2_matrix = tau2.as_matrix()
    tau2_matrix = np.delete(tau2_matrix,0,1)
    
    V2 = tau2_matrix.shape[0]
    G2 = tau2_matrix.shape[1]/4
    
    tau2_array = np.reshape(tau2_matrix,(V2, G2,4))
    
    comp = compSND(tau1_array,tau2_array)
    
    validate = np.min(comp,axis=1)/float(V2)
    
    np.savetxt(sys.stdout,validate)
if __name__ == "__main__":
    main(sys.argv[1:])