#!/usr/bin/env python3

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

from collections import defaultdict


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
    
    #import ipdb; ipdb.set_trace()
    tau_star_file = args.tau_star_file
    tau_file = args.tau_file

    tau_star = p.read_csv(tau_star_file, header=0, index_col=0)
    
    tau_star_var = defaultdict(dict)
    tau_star_matrix = tau_star.to_numpy()
    tau_star_matrix = np.delete(tau_star_matrix,0,1)
    tau_star_matrix[tau_star_matrix < 0.5] = 0.0
    tau_star_matrix[tau_star_matrix >= 0.5] = 1.0   
    
    idx = 0
    for gene, row in tau_star.iterrows():
        pos = row['Position']
        
        tau_star_var[gene][pos] = tau_star_matrix[idx,:]
        idx = idx + 1
    
    V = tau_star_matrix.shape[0]
    G = tau_star_matrix.shape[1]/4
    
    tau_star_array = np.reshape(tau_star_matrix,(V, G,4)) 
    
    tau = p.read_csv(tau_file, header=0, index_col=0)
    tau_matrix = tau.to_numpy()
    tau_matrix = np.delete(tau_matrix,0,1)
    
    PV = tau_matrix.shape[0]
    PG = tau_matrix.shape[1]/4
    
    #Want the intersection of the two sets

    idx = 0    
    intersect = 0
    tau_var = defaultdict(dict)
    
    for gene, row in tau.iterrows():
        pos = row['Position']
        
        tau_var[gene][pos] = tau_matrix[idx,:]
        if pos in tau_star_var[gene]:
            intersect += 1
        idx = idx + 1 
  
    tau_pred_intersect = np.zeros((intersect,G*4),dtype=int)
    tau_intersect = np.zeros((intersect,PG*4),dtype=int)
  
    idx = 0    
    intersect = 0
    
    for gene, row in tau.iterrows():
        pos = row['Position']
        
        tau_var[gene][pos] = tau_matrix[idx,:]
        if pos in tau_star_var[gene]:
            tau_pred_intersect[intersect,:] = tau_star_var[gene][pos]
            tau_intersect[intersect,:] = tau_var[gene][pos]
            intersect += 1
        idx = idx + 1 
        
  
    tau_star_pred = np.reshape(tau_pred_intersect,(intersect, G,4))
    tau_pred      = np.reshape(tau_intersect,(intersect, PG,4))
    
    #print "Intersection: " + str(intersect)
    
    comp = compSND(tau_star_pred,tau_pred)
    #print(comp)
    compV = comp/float(tau_star_pred.shape[0])

    accuracies = np.zeros(G)
    map = np.zeros(G,dtype=int)
    acctotal = 0.0
    ga = 0
    while (ga < G):
        (mr,mcol) = np.unravel_index(np.argmin(compV),compV.shape)
        curr_acc = np.min(compV)
        acctotal += curr_acc 
        compV[mr,:] = np.ones(PG)
        #compV[:,mcol] = np.ones(G)
        accuracies[mr] = curr_acc
        map[mr] = mcol
        ga += 1
        

    mean_acc = np.mean(accuracies)
    print(str(PG) + "," + str(G) + "," + str(mean_acc)) 

    #print(compV)
if __name__ == "__main__":
    main(sys.argv[1:])
