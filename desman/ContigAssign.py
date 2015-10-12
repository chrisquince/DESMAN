import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import argparse
import math
from operator import mul, div, eq, ne, add, ge, le, itemgetter

def get_sample_names(scg_freq):
    scgCols = scg_freq.columns.values.tolist()
    
    originalS = (len(scgCols) - 1)/4
    
    idx = range(1,originalS*4,4)
    
    sampleNames = [scgCols[i] for i in idx] 

    return sampleNames

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

    def __init__(self,cov,delta,n_run=None,max_iter=None,min_change=None):
        self.name = "KLAssign"

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
            self.min_change = 1.0e-3
        else:
            self.min_change = min_change
    
    def random_initialize(self):
        
        self.prng = np.random.RandomState()
            
        self.eta = self.prng.uniform(0, 1.0, (self.C, self.G))
        
            
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


def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("scg_freq_file", help="input core gene base frequencies")

    parser.add_argument("gamma_star_file", help="input MAP estimate frequencies")

    parser.add_argument("cov_file", help="mean contig coverages")

    parser.add_argument('-o','--output_stub', type=str, default="output",
        help=("string specifying output file stubs"))

    parser.add_argument('-g','--genomes', 
        help=("specify validation file of known genome composition."))

    args = parser.parse_args()

    import ipdb; ipdb.set_trace()

    scg_freq_file   = args.scg_freq_file
    gamma_star_file = args.gamma_star_file
    output_stub = args.output_stub
    cov_file = args.cov_file

    scg_freq    = p.read_csv(scg_freq_file, header=0, index_col=0)
    gamma_star = p.read_csv(gamma_star_file, header=0, index_col=0)
    cov = p.read_csv(cov_file, header=0, index_col=0)
    
    scg_freq_matrix = scg_freq.as_matrix()
    scg_freq_matrix = np.delete(scg_freq_matrix, 0, 1)
    
    scg_freqs = np.reshape(scg_freq_matrix, (scg_freq_matrix.shape[0],scg_freq_matrix.shape[1]/4,4))
    scg_totals = scg_freqs.sum(axis=2)
    
    sample_names = get_sample_names(scg_freq)
    scg_totals_df = p.DataFrame(scg_totals)
    scg_totals_df.columns = sample_names
     
    gamma_names = gamma_star.index.values
    
    
    scg_names = list(scg_totals_df.columns.values) 
    
    scg_names_new = [ rchop(x,'-A') for x in scg_names ]
    
    scg_totals_df.columns = scg_names_new
    
    scg_totals_df_gamma = scg_totals_df[gamma_names]
    
    scg_totals_gamma = scg_totals_df_gamma.as_matrix() 
    
    total_mean = np.mean(scg_totals_gamma,axis=0,dtype=np.float64)
    total_sd = np.std(scg_totals_gamma,axis=0,dtype=np.float64)
    
    mean_sd =  np.column_stack((total_mean,total_sd))
    mean_sd_df = p.DataFrame(mean_sd,index=gamma_names)
    mean_sd_df.columns = ['mean','sd']
    mean_sd_df.to_csv(output_stub+"mean_sd_df.csv")
    
    gamma_star_matrix = gamma_star.as_matrix()
    delta = np.multiply(gamma_star_matrix,total_mean[:,np.newaxis])
    
    cov = cov[gamma_names]
    cov_matrix = cov.as_matrix()
    klassign = KLAssign(cov_matrix,delta)
    klassign.factorize()
    
    contig_names = cov.index.tolist()
    etaD = np.rint(klassign.eta)
    etaD_df = p.DataFrame(etaD,index=contig_names)
    etaD_df.to_csv(output_stub+"etaD_df.csv")
    
    if args.genomes:
        genomes    = p.read_csv(args.genomes, header=0, index_col=0)
        genomes_M   = genomes.as_matrix()
        genomes_D = np.copy(genomes_M)
        
        genomes_D[genomes_D < 0.5] = 0.
        genomes_D[genomes_D >= 0.5] = 1.0
        
        NC = klassign.C
        NG = genomes_M.shape[1]
        acc = np.zeros((NG,klassign.G))
        for g in range(NG):
            for h in range(klassign.G):
                acc[g,h] = sum(genomes_D[:,g] == etaD[:,h])/float(NC)
                
        np.savetxt(sys.stdout, acc, fmt='%.4f')
        accuracies = np.zeros(NG)
        map = np.zeros(NG,dtype=int)
        acctotal = 0.0
        ga = 0;
        while (ga < NG):
            (mr,mcol) = np.unravel_index(np.argmax(acc),acc.shape)
            curr_acc = np.max(acc)
            acctotal += curr_acc 
            acc[mr,:] = np.zeros(klassign.G)
            acc[:,mcol] = np.zeros(NG)
            accuracies[mr] = curr_acc
            map[mr] = mcol
            ga += 1
        print str(acctotal) + "\n"
    
    
if __name__ == "__main__":
    main(sys.argv[1:])