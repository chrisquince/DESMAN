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
from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from numpy.random import RandomState

#user defined modules
import Desman_Utils as du

class Init_NMFT:
    """Initialises tau and gamma based on tensor non-negative matrix factorization""" 
   
    BASE_PRIOR = 1.0
   
    def __init__(self,snps, rank, randomState, n_run = 1, max_iter = 5000, min_change = 1.0e-5,alpha_constant=0.01):
        
        self.V = snps.shape[0] #number of variants
        self.S = snps.shape[1] # cos there are 4 bases
        self.G = rank;
        self.randomState = randomState 
        self.n_run = n_run;
        self.max_iter = max_iter;
        self.min_change = min_change;
        
        self.alpha = np.empty(self.G); self.alpha.fill(alpha_constant)
        self.alpha4 = np.empty(4); self.alpha4.fill(alpha_constant)
        
        variants = snps.astype(np.float) + self.BASE_PRIOR #prior on base frequency 
        
        #normalise frequencies at each site
        site_sum = variants.sum(axis = 2)
        freq = variants/site_sum[:,:,np.newaxis]
        
        self.N = self.V*4
        self.freq_matrix = np.zeros((self.N,self.S))
        
        for a in range(4):
            for v in range(self.V):
                self.freq_matrix[v + a*self.V,:] = freq[v,:,a] 
                
        self.tau = np.zeros((self.N,self.G))
        self.gamma = np.zeros((self.G,self.S))        
                
        
    def random_initialize(self):
        
        if self.G > 1:
            temp =  self.randomState.dirichlet(self.alpha, size=self.S)
        else:
            temp = np.ones((self.S,self.G))
        self.gamma = np.transpose(temp)
                
        for v in range(self.V):
            for g in range(self.G):
                tempvg = self.randomState.dirichlet(self.alpha4)
                for a in range(4):
                    self.tau[v + a*self.V,g] = tempvg[a] 
    
    def random_initialize_tau(self):
                
        for v in range(self.V):
            for g in range(self.G):
                tempvg = self.randomState.dirichlet(self.alpha4)
                for a in range(4):
                    self.tau[v + a*self.V,g] = tempvg[a] 
    
    def _adjustment(self):
        """Adjust small values to factors to avoid numerical underflow."""
        self.tau = np.maximum(self.tau, np.finfo(self.tau.dtype).eps)
        self.gamma = np.maximum(self.gamma, np.finfo(self.gamma.dtype).eps)
    
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
 
                if iter % 100 == 0:
                    logging.info('NTF Iter %d, div = %f'%(iter,div)) 

                iter += 1
                
    def factorize_gamma(self):
        for run in xrange(self.n_run):
           
            divl = 0.0
            div = self.div_objective()
            iter=0
            while iter < self.max_iter and math.fabs(divl - div) > self.min_change:
                self.div_update_gamma()
                #self._adjustment()
                divl = div
                div = self.div_objective()
 
                if iter % 100 == 0: 
                    print str(iter) + "," + str(div)

                iter += 1

    def factorize_tau(self):
        for run in xrange(self.n_run):
            self.random_initialize_tau()
            divl = 0.0
            div = self.div_objective()
            iter=0
            while iter < self.max_iter and math.fabs(divl - div) > self.min_change:
                self.div_update_tau()
                #self._adjustment()
                divl = div
                div = self.div_objective()
 
                if iter % 100 == 0:
                    logging.info('NTF Iter %d, div = %f'%(iter,div))

                iter += 1


    def div_objective(self):
        """Compute divergence of target matrix from its NMF estimate."""
        pa = du.elop(np.dot(self.tau, self.gamma))
        
        return (np.multiply(self.freq_matrix, np.log(du.elop(self.freq_matrix, pa, div))) - self.freq_matrix + pa).sum()
        
    def div_update(self):
        """Update basis and mixture matrix based on divergence multiplicative update rules."""
        H1 = np.tile(self.tau.sum(0)[:,np.newaxis],(1, self.S))
        
        if self.G > 1:
            self.gamma = np.multiply(self.gamma, du.elop(np.dot(self.tau.T, du.elop(self.freq_matrix, np.dot(self.tau, self.gamma), div)), H1, div))

            gamma_sum = self.gamma.sum(axis = 0)
            self.gamma = self.gamma/gamma_sum[np.newaxis,:]
        else:
            self.gamma = np.ones((self.G,self.S))

        tau1 = np.tile(self.gamma.sum(1)[np.newaxis,:],(self.N, 1))
        self.tau = np.multiply(
            self.tau, du.elop(np.dot(du.elop(self.freq_matrix, np.dot(self.tau, self.gamma), div), self.gamma.T), tau1, div)) 
            
        for v in range(self.V):
            for g in range(self.G):
                sumvg = 0.0
                for a in range(4):
                    sumvg += self.tau[v + a*self.V,g]
                    
                for a in range(4):
                    self.tau[v + a*self.V,g] = self.tau[v + a*self.V,g]/sumvg 

    def div_update_gamma(self):
        """Update basis and mixture matrix based on divergence multiplicative update rules."""
        H1 = np.tile(self.tau.sum(0)[:,np.newaxis],(1, self.S))
        
        self.gamma = np.multiply(self.gamma, du.elop(np.dot(self.tau.T, du.elop(self.freq_matrix, np.dot(self.tau, self.gamma), div)), H1, div))

        gamma_sum = self.gamma.sum(axis = 0)
        self.gamma = self.gamma/gamma_sum[np.newaxis,:]

    def div_update_tau(self):
        """Update basis and mixture matrix based on divergence multiplicative update rules."""
        tau1 = np.tile(self.gamma.sum(1)[np.newaxis,:],(self.N, 1))
        self.tau = np.multiply(
            self.tau, du.elop(np.dot(du.elop(self.freq_matrix, np.dot(self.tau, self.gamma), div), self.gamma.T), tau1, div)) 
            
        for v in range(self.V):
            for g in range(self.G):
                sumvg = 0.0
                for a in range(4):
                    sumvg += self.tau[v + a*self.V,g]
                    
                for a in range(4):
                    self.tau[v + a*self.V,g] = self.tau[v + a*self.V,g]/sumvg 



    def get_gamma(self):
        return np.transpose(self.gamma)

    def discretise_tau(self):
    
        #convert VX4 X G into VXGX4 shape
        discrete_tau = np.zeros((self.V*4,self.G))
        
        for v in range(self.V):
            for g in range(self.G):
                maxt = 0.0
                maxa = 0
                for a in range(4):
                    if(self.tau[v + a*self.V,g] > maxt):
                        maxt = self.tau[v + a*self.V,g] 
                        maxa = a 
                
                discrete_tau[v + self.V*maxa,g] = 1. 
        
        self.tau = discrete_tau
        
    def get_tau(self):
    
        #convert VX4 X G into VXGX4 shape
        ret_tau = np.zeros((self.V,self.G,4), dtype=np.int)
        
        for v in range(self.V):
            for g in range(self.G):
                maxt = 0.0
                maxa = 0
                for a in range(4):
                    if(self.tau[v + a*self.V,g] > maxt):
                        maxt = self.tau[v + a*self.V,g] 
                        maxa = a 
                ret_tau[v,g,maxa] = 1 
        
        return ret_tau