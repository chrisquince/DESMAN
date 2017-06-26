#!/usr/bin/env python

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
import logging

from operator import mul, div, eq, ne, add, ge, le, itemgetter
from itertools import izip
from itertools import compress
from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from numpy.random import RandomState
from scipy.stats import chi2
from collections import defaultdict

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
    new_pvalues = np.zeros(int(n))

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

def reject_outliers(data, m = 2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return s<m

class Variant_Filter():
    """Filters variant position based on simple binomial 
    or log ratio of binomial to mixtures of binomials"""

    def __init__(self,variants, randomState, optimise = True, threshold = 3.84, min_coverage = 5.0, 
                qvalue_cutoff = 0.1,max_iter = 100, min_p = 0.01, mCogFilter = 2.0,cogSampleFrac=0.95,Nthreshold=10):
        #first get array dimensions
        
        variants_matrix = variants.as_matrix()
        self.genes = list(variants.index)
        self.position = variants_matrix[:,0]
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
        self.mCogFilter = mCogFilter
        self.cogSampleFrac = cogSampleFrac
        self.max_iter = max_iter
        self.Nthreshold=Nthreshold
        self.eta = 0.96*np.identity((4)) + 0.01*np.ones((4,4))
        self.upperP = 1.0 - min_p
        #default select everything
        self.NS = self.V
        self.selected = np.ones((self.V), dtype=bool)    
        self.selected_indices  = list(np.where(self.selected))
        self.selected_indices = self.selected_indices[0].tolist()
    
        #denotes whether a random selection of significant positions has occured
        self.randomSelect = False
    
    def remove_outlier_cogs(self):
        totals = np.sum(self.freq,axis=1)
        
        uniquegenes = set(self.genes)
        
        gene_coverages = defaultdict(lambda: 0.)
        gene_freqs = defaultdict(lambda: 0.)
        idx = 0
        
        for gene in self.genes:
            gene_coverages[gene] += totals[idx]
            gene_freqs[gene] += 1.0
            idx = idx+1
        
        data = []
        for gene in uniquegenes:
            gene_coverages[gene] = gene_coverages[gene]/gene_freqs[gene]
            data.append(gene_coverages[gene])
        
        d = np.abs(data - np.median(data))
        mdev = np.median(d)
        s = d/mdev if mdev else 0.
        
        filterGenes = s > self.mCogFilter
    
        self.filteredGenes = [i for (i, v) in zip(list(uniquegenes), filterGenes) if v]
        
        idx = 0
        select = np.ones((self.V), dtype=bool)
        for gene in self.genes:
            if gene in self.filteredGenes:
                select[idx] = False
            
            idx+=1
    
        #stores filtered snp array
        self.snps_filter = self.snps_filter[select,:,:]
        
        self.V = self.snps_filter.shape[0] #number of variants
        
        #store base frequencies at each posn don't need sample info
        self.freq = (self.snps_filter).sum(axis=1) 
        self.ffreq = self.freq.astype(np.float)
        
        self.filtered = np.zeros((self.V), dtype=bool)

        self.NS = self.V
        self.genes = [i for (i, v) in zip(self.genes, select) if v]
        self.selected = np.ones((self.V), dtype=bool)    
        self.selected_indices  = list(np.where(self.selected))
        self.selected_indices = self.selected_indices[0].tolist()
        self.position = self.position[self.selected]
    
    
    def remove_outlier_cogs_sample(self):
        
        totals_samples = self.snps_filter.sum(axis=2)
        uniquegenes = list(set(self.genes))
        
        gene_sample_coverages = defaultdict(lambda: np.zeros(self.S))
        gene_sample_freqs = defaultdict(lambda: np.zeros(self.S))

        idx = 0
        
        for gene in self.genes:
            for s in range(self.S):
                gene_sample_coverages[gene][s] += totals_samples[idx][s]
                gene_sample_freqs[gene][s] += 1.0
            idx = idx+1
        
        nGenes = len(uniquegenes)
        
        g = 0
        geneSampleCovArray = np.zeros((nGenes,self.S))
        for gene in uniquegenes:
            for s in range(self.S):
                gene_sample_coverages[gene][s] = gene_sample_coverages[gene][s]/gene_sample_freqs[gene][s]
                geneSampleCovArray[g,s] = gene_sample_coverages[gene][s]
            g = g + 1    
            
        
        outlierGeneSample = np.zeros((nGenes,self.S),dtype=bool)
        for s in range(self.S):
            outlierGeneSample[:,s] = reject_outliers(geneSampleCovArray[:,s], m = self.mCogFilter)
    
        filterGenes = outlierGeneSample.sum(axis=1) < (self.S*self.cogSampleFrac)
        self.filteredGenes = [i for (i, v) in zip(uniquegenes, filterGenes) if v]
        
        idx = 0
        select = np.ones((self.V), dtype=bool)
        for gene in self.genes:
            if gene in self.filteredGenes:
                select[idx] = False
            
            idx+=1
    
        #stores filtered snp array
        self.snps_filter = self.snps_filter[select,:,:]
        
        self.V = self.snps_filter.shape[0] #number of variants
        
        #store base frequencies at each posn don't need sample info
        self.freq = (self.snps_filter).sum(axis=1) 
        self.ffreq = self.freq.astype(np.float)
        
        self.filtered = np.zeros((self.V), dtype=bool)

        self.NS = self.V
        self.genes = [i for (i, v) in zip(self.genes, select) if v]
        self.selected = np.ones((self.V), dtype=bool)
            
        self.selected_indices  = range(1,self.V)
        #list(np.where(select))
        #self.selected_indices = self.selected_indices[0].tolist()
        
        self.position = self.position[select]
    
    def selected_variants_todf(self,variants):

        contig_names = self.genes
    
        position = self.position
    
        snps_reshape = np.reshape(self.snps_filter,(self.NS,self.S*4))

        filtered_position = position[self.selected]
        selected_contig_names = [i for (i, v) in zip(contig_names, self.selected) if v]
    
        varCols = variants.columns.values.tolist()
        originalS = (len(varCols) - 1)/4
        sampleNames = list()
    
        j = 0
        for i in range(originalS):
            if j < self.S and i == self.sample_indices[j]:
            
                for a in range(4):
                    sampleNames.append(varCols[i*4 + 1 + a])
            
                j = j + 1
    
        snps_reshape_df = p.DataFrame(snps_reshape,index=selected_contig_names,columns=sampleNames)
    
        snps_reshape_df['Position'] = p.Series(filtered_position, index=snps_reshape_df.index)
        cols =  snps_reshape_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        snps_reshape_df = snps_reshape_df[cols]
        
        return snps_reshape_df
    
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
            ftemp[v,self.maxA[v]] = -1
        self.maxB = np.argmax(ftemp,axis=1)
        
        N = (self.freq.sum(axis=1)).astype(np.float64)
        n = (self.freq.max(axis=1)).astype(np.float64) #value of most abundant base
        m =  (ftemp.max(axis=1)).astype(np.float64)     #value of second most abundant
        e = N -n
        self.filtered = N < self.Nthreshold
        self.minV = np.zeros(self.V)
        self.minV[self.filtered == False] = m[self.filtered == False]/N[self.filtered == False] 
        
        p = np.zeros(self.V)
        MLL = np.zeros(self.V)
        lastSelect = 0
        Select = self.V
        p[self.filtered == False] = n[self.filtered == False]/N[self.filtered == False]
        p[p > self.upperP] = self.upperP
        while iter < self.max_iter and lastSelect != Select:
            #filter based on current error rate
        
            
            BLL = - (log(self.eta[self.maxA,:])*self.freq).sum(axis=1)
         
            for v in range(self.V):
                if  self.optimise:
                    res = minimize_scalar(mixNLL, bounds=(0.0, self.upperP), args = (self.eta,self.maxA[v],self.maxB[v],self.ffreq[v,:]), method='bounded')
                    p[v] = res.x
                MLL[v] = mixNLL(p[v],self.eta,self.maxA[v],self.maxB[v],self.ffreq[v,:])
            
            ratioNLL = 2.0*(BLL - MLL)
            
            self.filtered = np.logical_or(N < self.Nthreshold,ratioNLL < self.threshold)
            
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
            logging.info("Variant filter iter: " + str(iter) + " " + str(Select) + " " + str(self.eta))
            sys.stdout.flush()
            iter = iter + 1
            
        self.pvalue = 1.0 - chi2.cdf(ratioNLL,1)
        self.qvalue = benjamini_Hochberg(self.pvalue)
        self.ratioNLL = ratioNLL
        self.filtered = np.logical_or(N < self.Nthreshold,self.qvalue > self.qvalue_cutoff)
            
        self.snps_filter = self.snps_filter[self.filtered != True,:,:]
        self.selected_indices  = list(np.where(self.filtered != True))
        self.selected_indices = self.selected_indices[0].tolist()
        self.selected = self.filtered != True
        self.NS = self.snps_filter.shape[0]
        return self.snps_filter
        
    def select_Random(self, random_select):
        if(random_select < self.NS):
            self.randomSelect = True
        
            select = np.sort(self.randomState.choice(self.NS, random_select, replace=False))
            
            self.snps_filter_original = np.copy(self.snps_filter)
            self.snps_filter = self.snps_filter[select,:,:]
            
            self.NS = random_select
            
            self.selected_indices_original = np.copy(self.selected_indices)
            self.selected_indices = [self.selected_indices[i] for i in select]
        
            self.selected_original = np.copy(self.selected)
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
        help=("minimum coverage for sample to be included defaults 5.0"))
    
    parser.add_argument('-t','--outlier_thresh', type=float, default=2.0,
        help=("threshold for COG filtering on median coverage outlier defaults to 2.0"))

    parser.add_argument('-sf','--sample_frac', type=float, default=0.95,
        help=("fraction of samples with COG coverage exceeding median outlier for removal"))
    
    parser.add_argument('-o','--output_stub', type=str, default="output",
        help=("string specifying file stubs"))
        
    parser.add_argument('-p', '--optimiseP', action='store_true',help=("optimise proportions in likelihood ratio test default false"))
    
    parser.add_argument('-c', '--cog_filter', action='store_true',help=("whether to apply COG filtering default false"))
    
    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839"))
    
    args = parser.parse_args()
    
    variant_file = args.variant_file
    
    min_coverage = args.min_coverage
    output_stub = args.output_stub
    optimiseP = args.optimiseP
    random_seed = args.random_seed

    #start logging
    log_file_name = args.output_stub + 'log.txt'
    
    logging.basicConfig(
            filename=log_file_name,
            level=logging.INFO,
            filemode='w', # Overwrites old log file
            format='%(asctime)s:%(levelname)s:%(name)s:%(message)s'
            )

    logging.info("Results created at {0}".format(
            os.path.abspath(log_file_name)))


    #create new random state
    prng = RandomState(args.random_seed)
    
    max_qvalue = 1.0e-3
    if args.max_qvalue is not None:
        max_qvalue = args.max_qvalue
    
    filter_variants = 25.0 
    if args.filter_variants is not None:
        filter_variants = args.filter_variants
    
    min_variant_freq = 0.01
    if args.min_variant_freq is not None:
        min_variant_freq = args.min_variant_freq
        
    #read in snp variants
    #import ipdb; ipdb.set_trace()
    variants    = p.read_csv(variant_file, header=0, index_col=0)
    
    variant_Filter =  Variant_Filter(variants, randomState = prng, optimise = optimiseP, threshold = filter_variants, 
        min_coverage = min_coverage, qvalue_cutoff = max_qvalue, min_p = min_variant_freq, 
            mCogFilter = args.outlier_thresh,cogSampleFrac=args.sample_frac)
    
    if args.cog_filter:
        variant_Filter.remove_outlier_cogs_sample()
        filteredCogs_df = variant_Filter.selected_variants_todf(variants)
        filteredCogs_df.to_csv(output_stub+"cogf.csv")
    
    logging.info('Begun filtering variants with parameters: optimise probability = %s, lr threshold = %s, min. coverage = %s, q-value threshold = %s, min. variant frequency = %s' % (optimiseP, filter_variants, min_coverage, 
        max_qvalue, min_variant_freq))
    snps_filter = variant_Filter.get_filtered_VariantsLogRatio()
    logging.info("Completed variant filtering")       
    transition_matrix = variant_Filter.calc_Error_Matrix()
    
    contig_names = variant_Filter.genes
    
    position = variant_Filter.position
    
    snps_reshape_df = variant_Filter.selected_variants_todf(variants)
    
    logging.info('Output selected variants to %s',output_stub+"sel_var.csv") 
    snps_reshape_df.to_csv(output_stub+"sel_var.csv")
    
    logging.info('Output minor variant frequencies to %s',output_stub+"v_df.csv")
    p_df = p.DataFrame(variant_Filter.minV,index=contig_names)
    p_df = addPositions(p_df,position)
    p_df.to_csv(output_stub+"v_df.csv")

    logging.info('Output p-values to %s',output_stub+"p_df.csv")             
    p_df = p.DataFrame(variant_Filter.pvalue,index=contig_names)
    p_df = addPositions(p_df,position)   
    p_df.to_csv(output_stub+"p_df.csv")
    
    logging.info('Output q-values to %s',output_stub+"q_df.csv")
    q_df = p.DataFrame(variant_Filter.qvalue,index=contig_names)
    q_df = addPositions(q_df,position)  
    q_df.to_csv(output_stub+"q_df.csv")
        
    logging.info('Output log ratios to %s',output_stub+"r_df.csv")
    r_df = p.DataFrame(variant_Filter.ratioNLL,index=contig_names)
    r_df = addPositions(r_df,position)  
    r_df.to_csv(output_stub+"r_df.csv")
    
    logging.info('Output estimated error matrix to %s',output_stub+"tran_df.csv")
    etahat_df = p.DataFrame(transition_matrix)
    etahat_df.to_csv(output_stub+"tran_df.csv")
    
if __name__ == "__main__":
    main(sys.argv[1:])
