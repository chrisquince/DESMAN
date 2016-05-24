import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
from scipy.stats import norm
import argparse
import math
import Eta_Sampler2 as es
import sampletau

from operator import mul, div, eq, ne, add, ge, le, itemgetter
from numpy.random import RandomState

def intersect(a, b):
    return list(set(a) & set(b))


def get_sample_names(scg_freq):
    scgCols = scg_freq.columns.values.tolist()
    
    originalS = (len(scgCols) - 1)/4
    
    idx = range(1,originalS*4,4)
    
    sampleNames = [scgCols[i] for i in idx] 

    return sampleNames

def expand_sample_names(sample_names):
    expanded = []
    
    for name in sample_names:
        expanded.append(name + "-A")
        expanded.append(name + "-C")
        expanded.append(name + "-T")
        expanded.append(name + "-G")
    
    return expanded

def elop(Xt, Yt, op):
    X = np.copy(Xt)
    Y = np.copy(Yt)
    try:
        X[X == 0] = np.finfo(X.dtype).eps
        Y[Y == 0] = np.finfo(Y.dtype).eps
    except ValueError:
        return op(X, Y)
    return op(X, Y)
    
class KLAssign():

    def __init__(self,randomState,cov,delta,n_run=None,max_iter=None,min_change=None):
        self.name = "KLAssign"

        self.randomState = randomState
        self.cov = cov
        self.deltat = np.transpose(delta)

        self.C = self.cov.shape[0]
        self.G = self.deltat.shape[0]

        if n_run is None:
            self.n_run = 1
        else:
             self.n_run = n_run

        if max_iter is None:
            self.max_iter = 10000
        else:
             self.max_iter = max_iter

        if min_change is None:
            self.min_change = 1.0e-4
        else:
            self.min_change = min_change
    
    def random_initialize(self):
            
        self.eta = self.randomState.uniform(0, 1.0, (self.C, self.G))
        
            
    def factorize(self):
    
        for run in xrange(self.n_run):
            
            self.random_initialize()
        
            divl = 0.0
            div = self.div_objective()
            iter=0
            while iter < self.max_iter and math.fabs(divl - div) > self.min_change:
                self.div_update()
                self._adjustment()
                divl = div
                div = self.div_objective()
 
                print str(iter) + "," + str(div)

                iter += 1
    
    def _adjustment(self):
        """Adjust small values to factors to avoid numerical underflow."""
        self.eta = np.maximum(self.eta, np.finfo(self.eta.dtype).eps)
        
    
    def div_update(self):
        eta1 = np.tile(self.deltat.sum(1)[np.newaxis,:],(self.C, 1))
        self.eta = np.multiply(
            self.eta, elop(np.dot(elop(self.cov, np.dot(self.eta, self.deltat), div), self.deltat.T), eta1, div)) 
            
    def div_objective(self):
        """Compute divergence of target matrix from its NMF estimate."""
        cova = np.dot(self.eta, self.deltat)
        
        return (np.multiply(self.cov, np.log(elop(self.cov, cova, div))) - self.cov + cova).sum()

def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring

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

    parser.add_argument("scg_cov_file", help="input core gene coverages")

    parser.add_argument("gamma_star_file", help="input MAP estimate frequencies")

    parser.add_argument("cov_file", help="mean contig coverages")

    parser.add_argument("variants_file", help="variants called on contigs to be assigned")

    parser.add_argument("epsilon_file", help="predicted transition matrix")

    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839 applied after random filtering"))

    parser.add_argument('-o','--output_stub', type=str, default="output",
        help=("string specifying output file stubs"))

    parser.add_argument('-g','--genomes', 
        help=("specify validation file of known genome composition."))

    args = parser.parse_args()

    #import ipdb; ipdb.set_trace()

    prng = RandomState(args.random_seed)

    sampletau.initRNG()
    sampletau.setRNG(args.random_seed)
    output_stub = args.output_stub
    
    #read in data
    scg_cov    = p.read_csv(args.scg_cov_file, header=0, index_col=0)
    gamma_star = p.read_csv(args.gamma_star_file, header=0, index_col=0)
    cov = p.read_csv(args.cov_file, header=0, index_col=0)
    epsilon = p.read_csv(args.epsilon_file, header=0, index_col=0)
    epsilon_matrix = epsilon.as_matrix()
    
    variants = p.read_csv(args.variants_file, header=0, index_col=0)

    gamma_names = gamma_star.index.values
    scg_names = scg_cov.index.values

    intersect_names1 = sorted(intersect(gamma_names,scg_names))
    intersect_names = sorted(intersect(cov.columns.values,intersect_names1))
    
    scg_cov = scg_cov.reindex(intersect_names)
    gamma_star = gamma_star.reindex(intersect_names)
        
    total_mean = scg_cov['mean'].as_matrix()
    total_sd = scg_cov['sd'].as_matrix()
    
    #renormalise gamma matrix
    gamma_star_matrix = gamma_star.as_matrix()
    row_sums = gamma_star_matrix.sum(axis=1)
    gamma_star_matrix = gamma_star_matrix / row_sums[:, np.newaxis]

    delta = np.multiply(gamma_star_matrix,total_mean[:,np.newaxis])
    
    #reorder coverage matrix
    cov = cov[intersect_names]
    cov_matrix = cov.as_matrix()
    klassign = KLAssign(prng,cov_matrix,delta)
    klassign.factorize()
    
    #reorder variants matrix
    #var_sample_names = get_sample_names(variants)
    #var_sample_names_new = [ rchop(x,'-A') for x in var_sample_names ]
    #variants_cols = variants.columns.values
    
    expanded_names = expand_sample_names(intersect_names)
    
    variants = variants[expanded_names]
    
    #now apply Gaussian Gibbs sampler
    #import ipdb; ipdb.set_trace()
    etaD = np.rint(klassign.eta)
 
    etaSampler = es.Eta_Sampler2(prng,variants,cov,gamma_star_matrix,delta,total_sd,epsilon_matrix,etaD,max_var=20)
    etaSampler.update()
    contig_names = cov.index.tolist()
    
    etaD_df = p.DataFrame(etaD,index=contig_names)
    etaD_df.to_csv(output_stub+"etaD_df.csv")

    etaS_df = p.DataFrame(etaSampler.eta_star,index=contig_names)
    etaS_df.to_csv(output_stub+"etaS_df.csv")

    eta_df = p.DataFrame(klassign.eta,index=contig_names)
    eta_df.to_csv(output_stub+"eta_df.csv")   

    if args.genomes:
        genomes    = p.read_csv(args.genomes, header=0, index_col=0)
        genomes = genomes.loc[contig_names]
        genomes_M   = genomes.as_matrix()
        genomes_D = np.copy(genomes_M)
        
        genomes_D[genomes_D < 0.5] = 0.
        genomes_D[genomes_D >= 0.5] = 1.0
        
        etaD[etaD < 0.5] = 0.
        etaD[etaD >= 0.5] = 1.0
        
        (dtotal, dacc,dacc_array) = compGenes(etaD, genomes_D)
        (stotal, sacc,sacc_array) = compGenes(etaSampler.eta_star, genomes_D)
       
                
        np.savetxt(sys.stdout, dacc, fmt='%.4f')
        np.savetxt(sys.stdout, sacc, fmt='%.4f')
        
        print "KL: " + str(dtotal) + "\n"
        print "GS: " + str(stotal) + "\n"    
    
    
if __name__ == "__main__":
    main(sys.argv[1:])
