import glob
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
import re
from collections import defaultdict

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

def computeStrainReproducibility(gamma_file,tau_file,comp_files):

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
    if NC > 0:     
        mean_acc = np.mean(all_acc, axis=1)
    else:
        mean_acc = np.ones(G)
    #for g in range(G):
     #   print "%d,%f,%f" %(g,gamma_mean[g],mean_acc[g]) 

    return (gamma_mean,mean_acc)
def main(argv):

    #import ipdb; ipdb.set_trace()

    parser = argparse.ArgumentParser()
    
    parser.add_argument("input_stub", help="relative frequencies of haplotypes")
    
    parser.add_argument('-d','--delta_g', type=float, default=0.05,
        help=("minimum fractional reduction in PD default 0.05"))
    
    parser.add_argument('-m','--max_err', type=float, default=0.10,
        help=("maximum error valid strain"))
    
    parser.add_argument('-f','--min_freq', type=float, default=0.05,
        help=("minimum frequency valid strain"))
    
    args = parser.parse_args()
    
    deltaG = args.delta_g
    
    allPDG = defaultdict(defaultdict)
    
    #import ipdb; ipdb.set_trace()
    
    gValues = []
    
    for desman_dir in glob.glob(args.input_stub + "_*_0"):
       # print desman_dir
        
        if os.path.isdir(desman_dir):
            preg = re.compile('.*_(\d+)_0')
            m = preg.match(desman_dir)
            if m:
                G = int(m.group(1))
                gValues.append(G)
    
    gValues = sorted(gValues)
    NG      = len(gValues)
    sumPD   = np.zeros(NG)
    countPD = np.zeros(NG,dtype=np.int)
    
    gidx = 0
    for G in gValues:
        for fit_file in glob.glob(args.input_stub + "_" + str(G) + "_*/fit.txt"):
            p2 = re.compile(".*_" + str(G) + "_(\d+)")
            m2 = p2.match(fit_file)
            r = int(m2.group(1))
                    
            with open(fit_file, 'r') as f:
                fit_stat_line = f.readline().strip()
                (Ignore,GT,HT,LL,PD) = fit_stat_line.split(',')
                GT = int(GT)
                HT = int(HT)
                fPD = float(PD)
                #sometimes H smaller than G we are not going to use those replicates
                if HT == G:
                    allPDG[G][r] = fPD
                    sumPD[gidx]   += fPD
                    countPD[gidx] += 1                  
        gidx += 1
    
    meanPD = sumPD/countPD

    for gidx in range(2,NG):
        fractionalReduction = (meanPD[gidx - 1] - meanPD[gidx])/meanPD[gidx - 1] 

        if countPD[gidx] < 1 or fractionalReduction < deltaG:
            break
    
    newNG = gidx
    
    strainQuality = {}
    for gidx in range(newNG):
        G = gValues[gidx] 
        
        #get best replicate at this value of G
        bestr = min(allPDG[G], key=allPDG[G].get)
        gamma_file =  args.input_stub + "_" + str(G) + "_" + str(bestr) +"/Gamma_star.csv"
        tau_file = args.input_stub + "_" + str(G) + "_" + str(bestr) +"/Filtered_Tau_star.csv"
        comp_files = []
        
        for r, fPD in allPDG[G].iteritems():
            if r != bestr:
                comp_files.append(args.input_stub + "_" + str(G) + "_" + str(r) +"/Filtered_Tau_star.csv")
        (gamma_mean, mean_acc) = computeStrainReproducibility(gamma_file,tau_file,comp_files)
        
        hidx = 0
        NStrains = 0
        selected_err = []
        selected_strains = []
        for mean, acc in zip(gamma_mean, mean_acc):
            if mean > args.min_freq and acc < args.max_err:
                NStrains += 1
                selected_err.append(acc)
                selected_strains.append(hidx)
            hidx += 1
        if NStrains > 0:
            meanError = np.mean(selected_err)       
        else:
            #choose one strain in this event must abundant? most reproducible?
            bestidx = np.argmin(mean_acc)
            NStrains += 1
            selected_err.append(mean_acc[bestidx])
            selected_strains.append(bestidx)
            meanError = mean_acc[bestidx]
        
        strainQuality[G] = (NStrains,meanError,bestr,selected_strains,G)
    
    
    strainSorted = sorted(strainQuality, key=lambda k: (strainQuality[k][0], -strainQuality[k][1],-strainQuality[k][4]))
    
    
    if newNG > 0:
        bestG = strainSorted[-1]
        (NStrains,meanError,bestr,selected_strains,G) = strainQuality[bestG]
        
        tau_file = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Filtered_Tau_star.csv"
        
        print str(bestG)+ "," + str(NStrains) + "," + str(bestr) + "," + str(meanError) + "," + tau_file
        
        #Now remove low quality strains from gamma and tau files
        
        gamma_file = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Gamma_star.csv"
        
        gamma    = p.read_csv(gamma_file, header=0, index_col=0)
        
        gammaR = gamma.ix[:,selected_strains]
        
        gamma_fileR = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Gamma_starR.csv"
        
        gammaR.to_csv(gamma_fileR)
        
        
        
        gamma_file2 = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Gamma_mean.csv"
        
        gamma2    = p.read_csv(gamma_file2, header=0, index_col=0)
        
        gammaR2 = gamma2.ix[:,selected_strains]
        
        gamma_fileR2 = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Gamma_meanR.csv"
        
        gammaR2.to_csv(gamma_fileR2)
        
        
        
        tau = p.read_csv(tau_file, header=0, index_col=0)
        
        tau_select_idx = [0]
        for selected in selected_strains:
            offset = 1 + 4*selected
            for n in range(4):
                tau_select_idx.append(offset + n)
        
        tauR = tau.ix[:,tau_select_idx]
        
        tau_fileR = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Filtered_Tau_starR.csv"
        
        tauR.to_csv(tau_fileR)
        
        
        
        tau_mean_file = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Tau_Mean.csv"
        tau2 = p.read_csv(tau_mean_file, header=0, index_col=0)
        
        tauR2 = tau2.ix[:,tau_select_idx]
        
        tau_fileR2 = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Tau_MeanR.csv"
        
        tauR2.to_csv(tau_fileR2)
        
        tau_collated_file = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Collated_Tau_star.csv"
        
        if os.path.isfile(tau_collated_file):
            tauC = p.read_csv(tau_collated_file, header=0, index_col=0)
        
            tauCR = tauC.ix[:,tau_select_idx]
        
            tau_fileCR = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Collated_Tau_starR.csv"
        
            tauCR.to_csv(tau_fileCR)
            
            
            tau_collated_mean_file = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Collated_Tau_mean.csv"
        
            tauC2 = p.read_csv(tau_collated_mean_file, header=0, index_col=0)
        
            tauCR2 = tauC2.ix[:,tau_select_idx]
        
            tau_fileCR2 = args.input_stub + "_" + str(bestG) + "_" + str(bestr) +"/Collated_Tau_meanR.csv"
        
            tauCR2.to_csv(tau_fileCR2)
        
        #Gamma_star.csv  Tau_Mean.csv  
        #/Collated_Tau_star.csv, /Collated_Tau_mean.csv
        
           
if __name__ == "__main__":
    main(sys.argv[1:])
