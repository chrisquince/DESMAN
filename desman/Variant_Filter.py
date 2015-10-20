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
from numpy.random import RandomState
from scipy.stats import chi2

#class to perform simple variant filtering assuming fixed genome error rate

def log_factorial(x):
    """Returns the logarithm of x!
    Also accepts lists and NumPy arrays in place of x."""
    return gammaln(array(x)+1)

def mixNLL(p, eta, n, m, f):
    mix = p*eta[n,:] + (1-p)*eta[m,:]

    return np.dot(f,-np.log(mix))

def benjamini_Hochberg(pvalues):
    n = float(pvalues.shape[0])                                                                           
    new_pvalues = np.zeros(n)

    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]                                      
    values.sort()
    values.reverse()                                                                                  
    new_values = []
    for i, vals in enumerate(values):                                                                 
        rank = n - i
        pvalue, index = vals                                                                          
        new_values.append((n/rank) * pvalue)                                                          
    for i in xrange(0, int(n)-1):  
        if new_values[i] < new_values[i+1]:                                                           
            new_values[i+1] = new_values[i]                                                           
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i] 
        
    return new_pvalues
    
class Variant_Filter():
    """Filters variant position based on simple binomial 
    or log ratio of binomial to mixtures of binomials"""

    def __init__(self,variants, randomState, optimise = True, threshold = 3.84, min_coverage = 5.0, qvalue_cutoff = 0.1,max_iter = 100, min_p = 0.01):
        #first get array dimensions
        
        variants_matrix = variants.as_matrix()
        variants_matrix = np.delete(variants_matrix, 0, 1)
    
        snps = np.reshape(variants_matrix, (variants_matrix.shape[0],variants_matrix.shape[1]/4,4))
        vs_sum = snps.sum(axis=(2))
        vs_mean = np.mean(vs_sum,axis=0)
        #set random state
        self.randomState = randomState
        
        #stores samples filtered because min coverage too low
        
        self.sample_filter = vs_mean > min_coverage
        self.sample_indices  = list(np.where(self.sample_filter == True))
        self.sample_indices = self.sample_indices[0].tolist()
           
        #stores filtered snp array
        self.snps_filter = snps[:,self.sample_filter,:]
        
        self.V = self.snps_filter.shape[0] #number of variants
        
        self.S = self.snps_filter.shape[1] #number of sites following filtering
        
        #store base frequencies at each posn don't need sample info
        self.freq = (self.snps_filter).sum(axis=1) 
        self.ffreq = self.freq.astype(np.float)
        self.threshold = threshold
        self.qvalue_cutoff = qvalue_cutoff
        self.optimise = optimise
        self.filtered = np.zeros((self.V), dtype=bool)
        
        self.max_iter = max_iter
        
        self.eta = 0.96*np.identity((4)) + 0.01*np.ones((4,4))
        self.upperP = 1.0 - min_p
        #default select everything
        self.NS = self.V
        self.selected = np.ones((self.V), dtype=bool)    
        self.selected_indices  = list(np.where(self.selected))
        self.selected_indices = self.selected_indices[0].tolist()
    
    def get_filtered_Variants(self):
        iter = 0
        
        self.maxA = np.argmax(self.freq,axis=1)
        
        N = (self.freq.sum(axis=1)).astype(np.float64)
        n = (self.freq.max(axis=1)).astype(np.float64)
        e = N -n
        #self.maxN = 
        
        NLComb = -log_factorial(N) + log_factorial(self.freq).sum(axis=1)
        
        lastSelect = 0
        Select = self.V
        
        while iter < self.max_iter and lastSelect != Select:
            #filter based on current error rate
        
            NLL = NLComb - (log(self.eta[self.maxA,:])*self.freq).sum(axis=1)
            
            self.filtered = NLL < self.threshold
            
            ff = self.freq[self.filtered] 
            af = self.maxA[self.filtered]
            Nf = N[self.filtered]
            
            self.eta = 96*np.identity((4)) + np.ones((4,4))
            
            for v in range(Nf.shape[0]):
                self.eta[af[v],:] += ff[v,:]
                
            esums = self.eta.sum(axis = 1)
            self.eta = self.eta/esums[:,np.newaxis]
            
            lastSelect = Select
            Select = self.V - self.filtered.sum()
            #print str(iter) + " " + str(Select) + " " + str(self.eta)
            iter = iter + 1
        
        self.pvalue = 1.0 - chi2.cdf(NLL,1)
        self.qvalue = benjamini_Hochberg(self.pvalue)
        self.filtered = self.qvalue > self.qvalue_cutoff
        self.ratioNLL = NLL
        self.snps_filter = self.snps_filter[self.filtered != True,:,:]
        self.selected_indices  = list(np.where(self.filtered != True))
        self.selected_indices = self.selected_indices[0].tolist()
        self.selected = self.filtered != True
        self.NS = self.snps_filter.shape[0]
        return self.snps_filter
    
    def get_filtered_VariantsLogRatio(self):
        iter = 0
        
        #get most abundant base
        self.maxA = np.argmax(self.freq,axis=1)
        #get second most abundant base
        
        ftemp = np.copy(self.freq)
        for v in range(self.V):
            ftemp[v,self.maxA[v]] = 0
        self.maxB = np.argmax(ftemp,axis=1)
        
        N = (self.freq.sum(axis=1)).astype(np.float64)
        n = (self.freq.max(axis=1)).astype(np.float64) #value of most abundant base
        m =  (ftemp.max(axis=1)).astype(np.float64)     #value of second most abundant
        e = N -n
        
        p = np.zeros(self.V)
        MLL = np.zeros(self.V)
        lastSelect = 0
        Select = self.V
        
        p = n/N
        while iter < self.max_iter and lastSelect != Select:
            #filter based on current error rate
        
            #NLL = NLComb - (log(self.eta[self.maxA,:])*self.freq).sum(axis=1)
            
            BLL = - (log(self.eta[self.maxA,:])*self.freq).sum(axis=1)
         
            for v in range(self.V):
                if  self.optimise:
                    res = minimize_scalar(mixNLL, bounds=(0.0, self.upperP), args = (self.eta,self.maxA[v],self.maxB[v],self.ffreq[v,:]), method='bounded')
                    p[v] = res.x
                MLL[v] = mixNLL(p[v],self.eta,self.maxA[v],self.maxB[v],self.ffreq[v,:])
            
            ratioNLL = 2.0*(BLL - MLL)
            
            self.filtered = ratioNLL < self.threshold
            
            ff = self.freq[self.filtered] 
            af = self.maxA[self.filtered]
            Nf = N[self.filtered]
            
            self.eta = 96*np.identity((4)) + np.ones((4,4))
            
            for v in range(Nf.shape[0]):
                self.eta[af[v],:] += ff[v,:]
                
            esums = self.eta.sum(axis = 1)
            self.eta = self.eta/esums[:,np.newaxis]
            
            lastSelect = Select
            Select = self.V - self.filtered.sum()
            print str(iter) + " " + str(Select) + " " + str(self.eta)
            sys.stdout.flush()
            iter = iter + 1
            
        self.pvalue = 1.0 - chi2.cdf(ratioNLL,1)
        self.qvalue = benjamini_Hochberg(self.pvalue)
        self.ratioNLL = ratioNLL
        self.filtered = self.qvalue > self.qvalue_cutoff
            
        self.snps_filter = self.snps_filter[self.filtered != True,:,:]
        self.selected_indices  = list(np.where(self.filtered != True))
        self.selected_indices = self.selected_indices[0].tolist()
        self.selected = self.filtered != True
        self.NS = self.snps_filter.shape[0]
        return self.snps_filter
        
    def select_Random(self, random_select):
        if(random_select < self.NS):
            select = np.sort(self.randomState.choice(self.NS, random_select, replace=False))
        
            self.snps_filter = self.snps_filter[select,:,:]
            
            self.NS = random_select
            
            #all needs to be corrected to refer to original array.... 
            self.selected_indices = [self.selected_indices[i] for i in select]
        
            self.selected = np.zeros((self.V), dtype=bool)
            
            self.selected[self.selected_indices] = True
        
        return self.snps_filter
    
    def calc_Error_Matrix(self):
        self.tran_Matrix = np.ones((4,4)) #effective Laplace prior on transition matrix

        #select just sbv
        sbv_freq = self.freq[self.filtered]
        sbv_maxA = np.argmax(sbv_freq,axis=1)
        
        F = sbv_freq.shape[0]
        
        for i in range(F):
            self.tran_Matrix[sbv_maxA[i],] += sbv_freq[i]
            
        tSums = self.tran_Matrix.sum(axis=1)
        self.tran_Matrix = self.tran_Matrix / tSums[:,np.newaxis]

        return self.tran_Matrix

def addPositions(dataFrame,position):
    dataFrame['Position'] = position
    cols = dataFrame.columns.tolist()
    cols[0] = "value"
    dataFrame.columns = cols
    cols = cols[-1:] + cols[:-1]    
    dataFrame = dataFrame[cols] 
    return dataFrame
    
def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("variant_file", help="input SNP frequencies")
        
    parser.add_argument('-f','--filter_variants',nargs='?', const=3.84, type=float, 
        help=("binomial loge likelihood species p-value threshold for initial filtering as chi2"))
    
    parser.add_argument('-q','--max_qvalue',nargs='?', const=1.0e-3, type=float, 
        help=("specifies q value cut-off for variant defaults 1.0e-3"))
    
    parser.add_argument('-v','--min_variant_freq',nargs='?', const=0.01, type=float, 
        help=("specifies minimum variant frequency defaults 0.01"))
    
    parser.add_argument('-m','--min_coverage', type=float, default=5.0,
        help=("minimum coverage for sample to be included"))
    
    parser.add_argument('-o','--output_stub', type=str, default="output",
        help=("string specifying file stubs"))
        
    parser.add_argument('-p', '--optimiseP', action='store_true',help=("optimise proportions in likelihood ratio test"))
    
    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839"))
    
    args = parser.parse_args()
    
    variant_file = args.variant_file
    
    min_coverage = args.min_coverage
    output_stub = args.output_stub
    optimiseP = args.optimiseP
    random_seed = args.random_seed

    #create new random state
    prng = RandomState(args.random_seed)
    
    max_qvalue = 1.0e-3
    if args.max_qvalue is not None:
        max_qvalue = args.max_qvalue
    
    filter_variants = 3.84 
    if args.filter_variants is not None:
        filter_variants = args.filter_variants
    
    min_variant_freq = 0.01
    if args.min_variant_freq is not None:
        min_variant_freq = args.min_variant_freq
        
    #read in snp variants
    variants    = p.read_csv(variant_file, header=0, index_col=0)
    
    import ipdb; ipdb.set_trace()
    
    variant_Filter =  Variant_Filter(variants, randomState = prng, optimise = optimiseP, threshold = filter_variants, 
        min_coverage = min_coverage, qvalue_cutoff = max_qvalue, min_p = min_variant_freq)
    
    snps_filter = variant_Filter.get_filtered_VariantsLogRatio()
        
    transition_matrix = variant_Filter.calc_Error_Matrix()
    
    contig_names = variants.index.tolist()
    position = variants['Position']
    
    snps_reshape = np.reshape(variant_Filter.snps_filter,(variant_Filter.NS,variant_Filter.S*4))
    
    selected_Variants = variants[variant_Filter.selected]
    filtered_position = selected_Variants['Position']
    selected_contig_names = selected_Variants.index
    
    varCols = variants.columns.values.tolist()
    originalS = (len(varCols) - 1)/4
    sampleNames = list()
    
    j = 0
    for i in range(originalS):
        if j < variant_Filter.S and i == variant_Filter.sample_indices[j]:
            
            for a in range(4):
                sampleNames.append(varCols[i*4 + 1 + a])
            
            j = j + 1
    
    snps_reshape_df = p.DataFrame(snps_reshape,index=selected_contig_names,columns=sampleNames)
    
    snps_reshape_df['Position'] = filtered_position
    cols =  snps_reshape_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    snps_reshape_df = snps_reshape_df[cols]
    
    snps_reshape_df.to_csv(output_stub+"sel_var.csv")
                
    p_df = p.DataFrame(variant_Filter.pvalue,index=contig_names)
    p_df = addPositions(p_df,position)   
    p_df.to_csv(output_stub+"p_df.csv")
        
    q_df = p.DataFrame(variant_Filter.qvalue,index=contig_names)
    q_df = addPositions(q_df,position)  
    q_df.to_csv(output_stub+"q_df.csv")
        
    r_df = p.DataFrame(variant_Filter.ratioNLL,index=contig_names)
    r_df = addPositions(r_df,position)  
    r_df.to_csv(output_stub+"r_df.csv")
        
    etahat_df = p.DataFrame(transition_matrix)
    etahat_df.to_csv(output_stub+"tran_df.csv")
        
    with open(output_stub+"transition_matrix.pickle", 'w') as f:
        cPickle.dump(transition_matrix,f)
    
    with open(output_stub+"sel_var.pickle", 'w') as f:
        cPickle.dump(selected_Variants,f)
    
if __name__ == "__main__":
    main(sys.argv[1:])
