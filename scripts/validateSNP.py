

import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse
import pickle

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
        id0 = np.argmax(tau[v,0,:])
        for g in range(1,G):
            idg = np.argmax(tau[v,g,:])
            if(idg != id0):
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
                idg = np.argmax(tau1[v,g,:])
                idh = np.argmax(tau2[v,h,:])
                if(idg == idh):
                    overlap += 1 
                
            snd[g,h] = N - overlap
                
    return snd

def main(argv):

    #import ipdb; ipdb.set_trace()

    parser = argparse.ArgumentParser()
    parser.add_argument("tau_star_file", help="predicted variants")
        
    parser.add_argument("tau_file", help="known variants")
        
    args = parser.parse_args()
    
    tau_star_file = args.tau_star_file
    tau_file = args.tau_file

    tau_star = p.read_csv(tau_star_file, header=0, index_col=0)
    
    tau_star_matrix = tau_star.as_matrix()
    tau_star_matrix = np.delete(tau_star_matrix,0,1)
#    tau_star_matrix[tau_star_matrix < 0.5] = 0.0
#    tau_star_matrix[tau_star_matrix >= 0.5] = 1.0   

 
    V = tau_star_matrix.shape[0]
    G = tau_star_matrix.shape[1]/4
    
    tau_star_array = np.reshape(tau_star_matrix,(V, G,4)) 
    
    tau = p.read_csv(tau_file, header=0, index_col=0)
    tau_matrix = tau.as_matrix()
    tau_matrix = np.delete(tau_matrix,0,1)
    
    PV = tau_matrix.shape[0]
    PG = tau_matrix.shape[1]/4
    
    tau_array = np.reshape(tau_matrix,(PV, PG,4))
    
    #Want the intersection of the two sets
    taudict = {}
    idx = 0
    for index, row in tau.iterrows():
        taudict[(index,row['Position'])] = idx
        idx = idx + 1 
  
    var_true = np.zeros(V, dtype=bool)
    tau_pred_full = np.zeros((V,G,4),dtype=int)
    pindex = 0
    comp = 0
    for index, row in tau_star.iterrows():
        key = (index,int(row['Position']))
        if key in taudict: 
            comp = comp + 1
            var_true[pindex] = True
            tau_pred_full[pindex,:,:] =  tau_array[taudict[key],:,:]
        
        pindex = pindex + 1
    
    pred_true = variableTau(tau_star_array)
    
    correct_var = (pred_true == var_true)
    overlap = pred_true*var_true
    n_correct = correct_var.sum()
    p_correct = float(n_correct)/float(V)
    
    tau_star_pred = tau_star_array[overlap]
    tau_pred = tau_pred_full[overlap]
    
    
    comp = compSND(tau_star_pred,tau_pred)
    compV = comp/float(tau_star_pred.shape[0])

    print(compV)
if __name__ == "__main__":
    main(sys.argv[1:])
