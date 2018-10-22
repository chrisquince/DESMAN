#!/usr/bin/env python3
import argparse
import sys
import pandas as p
import numpy as np
from numpy.random import RandomState
from scipy.optimize import minimize_scalar
from numpy.random import RandomState
from scipy.stats import chi2
from collections import defaultdict
from scipy.special import gammaln

#class to perform simple variant filtering assuming fixed genome error rate
def div(x, y):
    if x.dtype == np.int64 and y.dtype == np.int64:
        return x // y
    return x / y

def log_factorial(x):
    """Returns the logarithm of x!
    Also accepts lists and NumPy arrays in place of x."""
    return gammaln(x + 1)

def mixNLL(p, eta, n, m, f):
    mix = p*eta[n,:] + (1-p)*eta[m,:]

    return np.dot(f,-np.log(mix))

def log_multinomial_pdf2(f,mix):
    return np.dot(f,-np.log(mix))
 

def log_multinomial_pdf(xs, ps):
    """Returns logarithm of multinomial pdf"""
    n = sum(xs)
    
    result = log_factorial(n) - sum(log_factorial(xs)) + sum(xs * np.log(ps))
    return result

def expand_sample_names(sample_names):
    expanded = []
    
    for name in sample_names:
        expanded.append(name + "-A")
        expanded.append(name + "-C")
        expanded.append(name + "-G")
        expanded.append(name + "-T")
    
    return expanded

class Variant_Val():
    """Filters variant position based on simple binomial 
    or log ratio of binomial to mixtures of binomials"""

    def __init__(self,variants, eta, G, tau_star_array, gamma_array, randomState, optimise = False, max_iter = 100, min_p = 0.01,sampleFilter=0.5,min_freq=10):
        self.optimise = optimise
        self.max_iter = max_iter
        self.min_p = min_p
        self.upperP = 1.0 - min_p
        self.min_freq = min_freq
        
        self.G = G
        self.tau = tau_star_array
        self.gamma = gamma_array
        
        row_sums = self.gamma.sum(axis=1)
        self.gamma = self.gamma / row_sums[:, np.newaxis]
        
        
        #first get array dimensions
        variants_matrix = variants.values
        self.genes = list(variants.index)
        self.position = variants_matrix[:,0]
        variants_matrix = np.delete(variants_matrix, 0, 1)
    
        self.snps = np.reshape(variants_matrix, (variants_matrix.shape[0],int(variants_matrix.shape[1] / 4),4))
        
        self.sumVS = np.sum(self.snps,axis=2)
        
        #self.snps = self.snps[:,np.sum(self.sumVS > 10,axis=0) > self.V*0.9]
        
        self.sumVS = np.sum(self.snps,axis=2)
        
        self.sampleFilter = sampleFilter
        self.S = self.snps.shape[1] #number of sites following filtering
        self.selectV = np.sum(self.sumVS > self.min_freq,axis=1) >= self.S*self.sampleFilter
        
        self.tau = self.tau[self.selectV,:,:]
        self.snps = self.snps[self.selectV,:]
        
        self.V = self.snps.shape[0] #number of variants
        
        self.sumVS = np.sum(self.snps,axis=2)
        self.position = self.position[self.selectV]
        self.genes = [x for x, y in zip(self.genes, self.selectV.tolist()) if y is True]
        
        self.snps1 = self.snps + 1
        self.fsnps = self.snps.astype(np.float)
        vs_sum = self.snps.sum(axis=(2))
        vs_mean = np.mean(vs_sum,axis=0)
        #set random state
        self.randomState = randomState
        self.eta = eta
        


    def calcMixtureLog(self):
        MLL = np.zeros(self.V)
        #get most abundant base
        self.maxA = np.argmax(self.snps,axis=2)
        self.maxAF = np.argmax(np.sum(self.snps,axis=1),axis=1)
        #get second most abundant base
        
        ftemp = np.copy(self.snps)
        for v in range(self.V):
            for s in range(self.S):
                ftemp[v,s,self.maxA[v,s]] = -1
        self.maxB = np.argmax(ftemp,axis=2)
        
        self.N = (self.snps.sum(axis=2)).astype(np.float64)
        self.n = (self.snps.max(axis=2)).astype(np.float64) #value of most abundant base
        self.m =  (ftemp.max(axis=2)).astype(np.float64)     #value of second most abundant
        self.e = self.N -self.n
        
        p = self.n/self.N
        p[p > self.upperP] = self.upperP
        q = self.m/self.N
        self.sumVS1 = np.sum(self.snps1,axis=2)
        
        self.nsnps = np.zeros((self.V,self.S,4))
        self.div = np.zeros(self.V)
        self.ndiv = np.zeros(self.V)
        
        for v in range(self.V):
            for s in range(self.S):
        
                #if self.sumVS[v,s] > 0:
                if self.sumVS[v,s] > self.min_freq:
                    fracd = self.snps1[v,s,:]/self.sumVS1[v,s]
                
                    self.div[v] += np.sum(self.probVS[v,s,:]*np.log(self.probVS[v,s,:]/fracd))
        
                    self.ndiv[v] += 1
        
                if  self.optimise:
                    res = minimize_scalar(mixNLL, bounds=(0.0, self.upperP), args = (self.eta,self.maxA[v,s],self.maxB[v,s],self.fsnps[v,s,:]), method='bounded')
                    p[v,s] = res.x

                    MLL[v] += mixNLL(p[v,s],self.eta,self.maxA[v,s],self.maxB[v,s],self.fsnps[v,s,:])
                else:
                    if self.sumVS[v,s] > self.min_freq:
                        frac1 = self.snps1[v,s,:]/self.sumVS1[v,s]
                        
                        self.nsnps[v,s,:] = self.snps[v,s,:]/self.sumVS[v,s]
                        #self.nsnps[v,s,:] = self.snps1[v,s,:]/self.sumVS1[v,s]
                        mix = np.dot(frac1,self.eta)
                    #p[v,s]*self.eta[self.maxA[v,s],:] + q[v,s]*self.eta[self.maxB[v,s],:]
                    
                        MLL[v] += np.dot(self.fsnps[v,s,:],-np.log(mix))
                        
        return MLL

    def calcRankTable(self):
    
        self.rsnps = np.zeros((self.V,self.S,4))
        
        for v in range(self.V):
            for s in range(self.S):
                self.rsnps[v,s,:] = np.sort(self.nsnps[v,s,:])
    
    def calcRankTable1(self):
        
        self.rsnps1 = np.zeros((self.V,self.S))
        
        for v in range(self.V):
            self.rsnps1[v,:] = 1.0 - self.nsnps[v,:,self.maxAF[v]]
                
    
    def calcModeLog(self):
        
        self.probVS = np.einsum('ijk,lj,km->ilm',self.tau,self.gamma,self.eta)
        sumVS = np.sum(self.snps,axis=2)
        logProbV = np.zeros(self.V)
        #loop each variant
        for v in range(self.V):
            #loop each site
            for s in range(self.S):
                if sumVS[v,s] > 0:                   
                    logProbV[v] += np.dot(self.snps[v,s,:],-np.log(self.probVS[v,s,:]))

        return logProbV

    
    def calcAssArray(self):
    
        self.assArray = np.zeros(self.V,dtype=np.int)
    
        self.tauT = np.zeros((self.V,self.G),dtype=np.int)
        for v in range(self.V):
            for g in range(self.G):
                self.tauT[v,g] = np.argmax(self.tau[v,g,:])
    
    
        self.mapAss = {}
        
        
        
        for v in range(self.V):
            test = (self.tauT[v,:] == self.maxAF[v]) == False
            
            testI = test.astype(int)
            testBit = np.packbits(testI)
            assert 8 - self.G > 0
            testBitR = np.right_shift(testBit, 8 - self.G)
        
            if np.sum(test) < self.G:
                self.assArray[v] = testBitR[0]
    
def main(argv):
    parser = argparse.ArgumentParser()
    
    parser.add_argument("variant_file", help="input SNP frequencies")

    parser.add_argument("tau_file", help="haplotype assignments")

    parser.add_argument("gamma_file", help="haplotype frequencies")

    parser.add_argument("eta_file", help="error transition matrix")

    parser.add_argument("outStub", help="output file stub")

    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839"))

    parser.add_argument('-m','--min_freq', type=int, default=5,
        help=("minimum frequency for sample to be used in SNV stats"))
    
    parser.add_argument('-f','--sample_frac', type=float, default=0.5,
        help=("fraction of samples with freq above minimum for SNV to be used at all"))

    parser.add_argument('-d','--min_dev', type=float, default=10,
        help=("number of samples for us to trust SNV stats"))

    args = parser.parse_args()
    
    #create new random state
    prng = RandomState(args.random_seed)
    
    
    variant_file = args.variant_file
    
    tau_star = p.read_csv(args.tau_file, header=0, index_col=0)
    
    tau_star_matrix = tau_star.as_matrix()
    tau_star_matrix = np.delete(tau_star_matrix,0,1)
    tau_star_matrix[tau_star_matrix < 0.5] = 0.0
    tau_star_matrix[tau_star_matrix >= 0.5] = 1.0   
        
    V = tau_star_matrix.shape[0]
    G = int(tau_star_matrix.shape[1]/4)
    
    tau_star_array = np.reshape(tau_star_matrix,(V, G,4)) 

    #import ipdb; ipdb.set_trace()
    
    variants    = p.read_csv(variant_file, header=0, index_col=0)
    
    eta_df = p.read_csv(args.eta_file, header=0, index_col=0)
    eta = eta_df.values
    
    gamma = p.read_csv(args.gamma_file, header=0, index_col=0)
    
    gamma_names = gamma.index.values
    expanded_names = expand_sample_names(gamma_names)
    variants_intersect = variants[['Position'] + expanded_names]
    
    gamma_array = gamma.values 
    
    variantVal = Variant_Val(variants_intersect,eta,G,tau_star_array,gamma_array,randomState = prng,sampleFilter=args.sample_frac,min_freq=args.min_freq)
    
    logModel = variantVal.calcModeLog()
    
    logMixture = variantVal.calcMixtureLog()
    
    logDiff = logMixture - logModel
    
    variantVal.calcRankTable1()

    df = p.DataFrame(variantVal.rsnps1,index=variantVal.position)
    df.to_csv(args.outStub + "_rsnps1.csv")
    
    variantVal.calcAssArray()
    
    with open(args.outStub + "_dev.csv", 'w') as out:
        for v in range(variantVal.V):
            maxD = max(-logDiff[v],0)
            if variantVal.ndiv[v] > 0:
                out.write(str(v) + "," + variantVal.genes[v] + "," + str(variantVal.position[v]) + "," + str(variantVal.ndiv[v]) + "," + str(variantVal.assArray[v]) +"," + str(maxD/variantVal.ndiv[v]) + "," + str(variantVal.div[v]/variantVal.ndiv[v]) + "\n")
            else:
                out.write(str(v) +  "," + variantVal.genes[v] + "," + str(variantVal.position[v]) + ","  + str(variantVal.ndiv[v]) + "," + str(variantVal.assArray[v]) +"," + "NA" + "," + "NA\n")
if __name__ == "__main__":
    main(sys.argv[1:])
