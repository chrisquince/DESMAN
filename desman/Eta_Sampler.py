import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
from scipy.stats import norm
import math
import sampletau
from scipy.special import gammaln
from numpy import array, log, exp
import Init_NMFT as inmft
import logging

MIN_DELTA = 1.0e-10
ETA_PENALTY=-1.0e3

def log_factorial(x):
    """Returns the logarithm of x!
    Also accepts lists and NumPy arrays in place of x."""
    return gammaln(array(x)+1)


def log_multinomial_pdf(xs, ps):
    """Returns logarithm of multinomial pdf"""
    n = sum(xs)
    xs, ps = array(xs), array(ps)
    result = log_factorial(n) - sum(log_factorial(xs)) + sum(xs * log(ps))
    return result
    
def log_Poisson(cov,cov_exp):
    
    cov_exp[cov_exp < MIN_DELTA] = MIN_DELTA    
    res = -gammaln(array(cov) + 1.) - cov_exp + cov*np.log(cov_exp)
        
    return res;

class Eta_Sampler():
    
    def __init__(self,randomState,variants,covs,gamma,delta,cov_sd,epsilon,init_eta,max_iter=None,max_eta=2,eta_scale=0.01,max_var=None):
    
        #calc G
        self.randomState = randomState
        self.delta = np.transpose(delta)
        self.cov_sd = np.transpose(cov_sd)
        self.gamma = np.copy(gamma,order='C')
        self.cov = covs.as_matrix()
        self.epsilon = np.copy(epsilon,order='C')
        
        self.G = self.gamma.shape[1]
        self.C = self.cov.shape[0]
        self.S = self.gamma.shape[0]
        
        
        self.ll = 0.0
        self.gene_ll = np.zeros(self.C)
        self.gene_llstar = np.zeros(self.C)
        
        self.genes = covs.index.tolist()
        
        self.gene_variants = {}
        self.gene_variants_full = {}
        self.gene_V = {}
        self.gene_tau = {}
        self.gene_map = {}
        
        c = 0
        for gene in self.genes:
            
            try:
                if variants is None:
                    raise KeyError()
                
                gene_variants = variants.loc[gene]
                
                #print gene + " " + str(gene_variants.shape[0])
                if len(gene_variants.shape) == 1:
                    gene_variants = gene_variants.to_frame()
                    gene_variants = gene_variants.transpose()
            
                gene_variants_matrix = gene_variants.as_matrix()
                
                gene_snps = np.reshape(gene_variants_matrix, (gene_variants_matrix.shape[0],gene_variants_matrix.shape[1]/4,4))
                NV = gene_snps.shape[0]
                if max_var is not None and NV > max_var:
                    select = np.sort(self.randomState.choice(NV,max_var, replace=False))
                    self.gene_variants_full[gene] = gene_snps.astype(int,order='C')
                    gene_snps = gene_snps[select,:,:]
                    
                self.gene_variants[gene] = gene_snps.astype(int,order='C')
                self.gene_V[gene] = gene_snps.shape[0]
                
                self.gene_tau[gene] = np.zeros((self.gene_V[gene],self.G,4), dtype=np.int,order='C')
                
                pass
            except KeyError:
                self.gene_variants[gene] = None
                self.gene_V[gene] = 0
                self.gene_tau[gene] = None
                
                pass
                        
            self.gene_map[gene] = c 
            
            c = c +1
    
        if max_iter is None:
            self.max_iter = 20
        else:
            self.max_iter = max_iter
            
        self.max_eta = max_eta 
        self.eta = np.copy(init_eta)
        self.eta[self.eta > self.max_eta - 1.0] = self.max_eta - 1.0

        self.eta_star = np.copy(init_eta)
        self.eta_store = np.zeros((self.max_iter,self.C,self.G))
        
        for gene in self.genes:
            c = self.gene_map[gene]
            if self.gene_V[gene] > 0:
                init_NMFT = inmft.Init_NMFT(self.gene_variants[gene],self.G,self.randomState)
                gammaR = self.maskGamma(self.gamma,self.eta[c,:])
                init_NMFT.gamma = np.transpose(gammaR)
                init_NMFT.factorize_tau()
                self.gene_tau[gene] = np.copy(init_NMFT.get_tau(),order='C')
                self.sampleTauC(self.gene_tau[gene],self.gene_variants[gene],self.eta[c,:])    
    
    

        self.eta_scale = eta_scale
        self.eta_log_prior = np.zeros(self.max_eta)
        for s in range(self.max_eta):
            self.eta_log_prior[s] = s*np.log(self.eta_scale)
        
        tempnorm = np.sum(np.exp(self.eta_log_prior))
        self.eta_log_prior -= np.log(tempnorm)

    def maskGamma(self,gamma,eta):
        gammaR = np.copy(gamma)
        
        for g in range(self.G):
            if eta[g] == 0:
                gammaR[:,g] = 0.0
        
        row_sums = gammaR.sum(axis=1)
        gammaR = gammaR / row_sums[:, np.newaxis]
    
        return gammaR
        
        
    def logLikelihoodGene(self,variants,cTau,cGamma,eta,cEpsilon):
        """Computes data log likelihood given parameter states"""
        logLL = 0.0
        V = variants.shape[0]
        if eta.sum() > 0:
            gammaR = self.maskGamma(cGamma,eta)
        
            probVS = np.einsum('ijk,lj,km->ilm',cTau,gammaR,cEpsilon)
        
            #loop each variant

            for v in range(V):
            #loop each site
                for s in range(self.S):                    
                    logLL += log_multinomial_pdf(variants[v,s,:], probVS[v,s,:])
        else:
            logLL = float(V)*ETA_PENALTY
        return logLL
    
    def logLikelihood(self):
    
        cov_expminus = np.dot(self.eta,self.delta)
        
        sumLogLL = 0.
       
        for gene in self.genes:
            c = self.gene_map[gene]
            self.gene_ll[c] = 0.0
            for g in xrange(self.G):
                self.gene_ll[c] += self.eta_log_prior[self.eta[c,g]]
            
            self.gene_ll[c] += log_Poisson(self.cov[c,:],cov_expminus[c,:]).sum()
            
            if self.gene_V[gene] > 0:
                self.gene_ll[c] += self.logLikelihoodGene(self.gene_variants[gene],self.gene_tau[gene],self.gamma,self.eta[c,:],self.epsilon)    
            
            #print str(c) + "," + str(self.gene_ll[c])  
            
            sumLogLL += self.gene_ll[c]
        return sumLogLL
    
    def computeVarLLContrib(self,eta,tau,variants):
        #then with strain removed
        gammaT = self.maskGamma(self.gamma,eta)
                            #definitely update tau for added gene
        newTau = np.copy(tau)
        nchange = self.sampleTauC(newTau,variants,eta)
                        
        siteVariantsProb = np.einsum('ijk,lj,km->ilm',newTau,gammaT,self.epsilon)
        st0 = np.log(siteVariantsProb)*variants    
        logprob0 = st0.sum()
    
        return(logprob0, newTau)
    
    def update(self): #perform max_iter Gibbs updates
        iter = 0
        self.ll = self.logLikelihood()
        for c in xrange(self.C):
            self.eta_star[c,:] = np.copy(self.eta[c,:])
            self.gene_llstar[c] = self.gene_ll[c]
        
        while (iter < self.max_iter):
        
            for gene in self.genes:
                c = self.gene_map[gene]
                V = self.gene_V[gene]
                
                
                for g in range(self.G):
                    
                    #first add coverages
                    state_logprob = np.copy(self.eta_log_prior)
                
                    tempEta = np.copy(self.eta[c,:])
                    tempEta[g] = 0

                    logVar0 = 0.

                    if V > 0:
                        variants = self.gene_variants[gene]
                        if(tempEta.sum() > 0):
                            (logVar0, newTau0) = self.computeVarLLContrib(tempEta,self.gene_tau[gene],variants)
                        else:
                            logVar0 = -1.0e20;
                        
                    
                    cov_expminus = np.dot(tempEta,self.delta)
                    #temp = norm.logpdf(self.cov[c,:],cov_expminus, self.cov_sd)
                    temp = log_Poisson(self.cov[c,:],cov_expminus)
                    state_logprob[0] += temp.sum() + logVar0
                    
                    tempEta[g] = 1.
                    
                    if V > 0:
                        variants = self.gene_variants[gene]
                        (logVar1, newTau1) = self.computeVarLLContrib(tempEta,self.gene_tau[gene],variants)
                    else:
                        logVar1 = 0.0
                    
                    for s in xrange(1,self.max_eta):
                        cov_expminus_S = cov_expminus + s*self.delta[g,:]
                        #temp = norm.logpdf(self.cov[c,:],cov_expminus, self.cov_sd)
                        temp = log_Poisson(self.cov[c,:],cov_expminus_S)
                        state_logprob[s] += temp.sum() + logVar1
                
                    eta_sample = self.sampleLogProb(state_logprob)
                    self.eta[c,g] = eta_sample
                    
                    if V > 0:
                        if eta_sample == 0:
                            self.gene_tau[gene] = newTau0
                        else:
                            self.gene_tau[gene] = newTau1
                        
            self.ll = self.logLikelihood()
   
            logging.info('Gibbs Iter %d, nll = %f'%(iter,self.ll))
            
            self.storeStarState(iter)
            self.eta_store[iter,]=np.copy(self.eta)

            iter = iter + 1
    
    def update2(self): #perform max_iter Gibbs updates
        iter = 0
        self.ll = self.logLikelihood()
        for c in xrange(self.C):
            self.eta_star[c,:] = np.copy(self.eta[c,:])
            self.gene_llstar[c] = self.gene_ll[c]
        
        while (iter < self.max_iter):

                
            for gene in self.genes:
                c = self.gene_map[gene]
                V = self.gene_V[gene]
                variants = self.gene_variants[gene]
                
                for g in range(self.G):
                    
                    #first add coverages
                    state_logprob = np.copy(self.eta_log_prior)
                    
                    tempEta = np.copy(self.eta[c,:])
                    newTaus = np.zeros((self.max_eta,V,self.G,4), dtype=np.int,order='C')
                    for s in xrange(0,self.max_eta):
                        tempEta[g] = s
                        
                        logVar = 0.0
                        
                        if V > 0:
                            newTaus[s,:] = self.gene_tau[gene]
                            if tempEta.sum() > 0:
                                (logVar, newTaus[s,:]) = self.computeVarLLContrib(tempEta,self.gene_tau[gene],variants)
                            else:
                                logVar = -1.0e20;
                        
                        cov_expminus = np.dot(tempEta,self.delta)
                        
                        temp = log_Poisson(self.cov[c,:],cov_expminus)
                            
                        state_logprob[s] += temp.sum() + logVar
                
                    eta_sample = self.sampleLogProb(state_logprob)
                    
                    self.eta[c,g] = eta_sample
                    
                    if V > 0:
                        self.gene_tau[gene] = newTaus[s]
                        
            self.ll = self.logLikelihood()
            print "Iter = %d: ll = %f\n" %(iter,self.ll)
            
            self.storeStarState(iter)
            self.eta_store[iter,]=np.copy(self.eta)

            
            #if (iter % 10 == 0):    
             #   logging.info('Gibbs Iter %d, no. changed = %d, nll = %f'%(iter,nchange,self.ll))
            
            iter = iter + 1
    

    def baseProbabilityGivenTau(self,tauState,gamma,epsilon):
        #returns 4 dim vector of each base probability A,C,G,T
        #given gamma and defined assignment tauState
        #loop bases

        return np.einsum('jk,lj,km->lm',tauState,gamma,epsilon)

    def sampleLogProb(self,adLogProbS):
        dP = np.exp(adLogProbS - np.max(adLogProbS))
        dP = dP/np.sum(dP,axis=0)
        return np.flatnonzero(self.randomState.multinomial(1,dP,1))[0]

    def sampleTauC(self,tau,variants,eta,gamma=None,epsilon=None):

        if gamma is None:
            gamma = self.gamma
        
        if epsilon is None:
            epsilon = self.epsilon
        
        V = tau.shape[0]
        
        #rework gamma matrix
        gammaR = self.maskGamma(gamma,eta)
        nchange = sampletau.sample_tau(tau, gammaR, epsilon, variants)

        return nchange
    
    def restoreFullVariants(self):
        for gene_full in self.gene_variants_full:
            self.gene_variants[gene] = self.gene_variants_full[gene]

            self.gene_V[gene] = self.gene_variants[gene].shape[0]
                
            self.gene_tau[gene] = np.zeros((self.gene_V[gene],self.G,4), dtype=np.int,order='C')
    
    def tauLikelihoodGene(self,variants,cTau,cGamma,eta,cEpsilon):
        """Computes data log likelihood given parameter states"""

        V = variants.shape[0]
        logLL = np.zeros(V)
        if eta.sum() > 0:
            gammaR = self.maskGamma(cGamma,eta)
        
            probVS = np.einsum('ijk,lj,km->ilm',cTau,gammaR,cEpsilon)
        
            #loop each variant

            for v in range(V):
                for s in range(self.S):                    
                    logLL[v] += log_multinomial_pdf(variants[v,s,:], probVS[v,s,:])
        
        return logLL
    
    def calcTauStar(self,eta,gamma=None,epsilon=None):
        if gamma is None:
            gamma = self.gamma
        
        if epsilon is None:
            epsilon = self.epsilon

        iter = 0
        
        self.gene_tau_star = {}
        self.gene_ll_tau_star = {}
        
        for gene in self.genes:
            V = self.gene_V[gene]
            self.gene_ll_tau_star[gene] = np.zeros(V)
            self.gene_ll_tau_star[gene].fill(np.finfo(np.float).min)
            self.gene_tau_star[gene] = np.zeros((V,self.G,4), dtype=np.int,order='C')
         
        while (iter < self.max_iter):
            lltausum = 0.0
            for gene in self.genes:
                c = self.gene_map[gene]
                V = self.gene_V[gene]
                etaSum = eta[c,:].sum()
                
                if V > 0:
                    if etaSum > 0:
                        nchange = self.sampleTauC(self.gene_tau[gene],self.gene_variants[gene],eta[c,:],gamma,epsilon)
                        tauLL = self.tauLikelihoodGene(self.gene_variants[gene],self.gene_tau[gene],gamma,eta[c,:],epsilon)
                        
                        for v in range(V):
                            if tauLL[v] > self.gene_ll_tau_star[gene][v]:
                                self.gene_ll_tau_star[gene][v] = tauLL[v]
                                self.gene_tau_star[gene][v,:] = np.copy(self.gene_tau[gene][v,:],order='C')
                        
                        lltausum += self.gene_ll_tau_star[gene].sum()
                    else:
                        nchange = -1
                    logging.info('Tau star Iter %d, nll = %f'%(iter,lltausum))
                    #print "c = " + str(c) + ", change = " + str(nchange) 
            
            
            iter = iter + 1
            #print "Iter = " + str(iter) + ", ll = " + str(lltausum) 
            
    def sampleTau(self,tau,variants,eta,gamma=None,epsilon=None):

        if gamma is None:
            gamma = self.gamma
        
        if epsilon is None:
            epsilon = self.epsilon
        
        V = tau.shape[0]
        
        #rework gamma matrix
        gammaR = self.maskGamma(gamma,eta)
        #nchange = sampletau.sample_tau(tau, gammaR, epsilon, variants)
        nchange = 0
        for v in range(V):
            #calculate probability of assignment of each genome to 1 of 4 bases
            for g in range(self.G):
                if eta[g] > 0.0:
                    propTau = np.zeros((4,self.G,4),dtype=np.int)    
                    stateLogProb = np.zeros(4)
                
                    for a in range(4):
                        propTau[a,:,:] = tau[v,:,:]
                        propTau[a,g,:] = np.zeros(4,dtype=np.int)
                        propTau[a,g,a] = 1
                        
                        siteProb  = self.baseProbabilityGivenTau(propTau[a,:,:],gammaR,epsilon)
                        st1 = np.log(siteProb)*variants[v,:,:]
                        stateLogProb[a] = st1.sum()
                
                        s = self.sampleLogProb(stateLogProb)
                
                    if not np.array_equal(propTau[s,:,:],tau[v,:,:]):
                        tau[v,:,:] = propTau[s,:,:]
                        nchange+=1
                    
        return nchange
    
    def getTauStar(self,variants):
        C = len(self.genes)
        Varray = np.zeros(C,dtype=np.int,order='C')  
        Vtotal = 0
        #vstart = {};
        
        for gene in self.genes:
            c = self.gene_map[gene]
            Varray[c] = self.gene_V[gene]    
            Vtotal += Varray[c]
        
        Vcum_array = np.zeros(C,dtype=np.int,order='C')
        for c in range(1,C):
            Vcum_array[c] = Vcum_array[c - 1] + Varray[c - 1]
        
        tauStar =  np.zeros((Vtotal,self.G,4), dtype=np.int,order='C')   
    
        for gene in self.genes:
            c = self.gene_map[gene]
            V = self.gene_V[gene]
            start = Vcum_array[c]
            end = start + V
            tauStar[start:end,:] = self.gene_tau_star[gene]
    
        positions = np.zeros(Vtotal,dtype=np.int,order='C')
        contig_index = ["" for x in range(Vtotal)]
        
        for gene in self.genes:
            c = self.gene_map[gene]
            V = self.gene_V[gene]
        
        
            try:
                gene_variants = variants.loc[gene]
                
                #print gene + " " + str(gene_variants.shape[0])
                if len(gene_variants.shape) == 1:
                    gene_variants = gene_variants.to_frame()
                    gene_variants = gene_variants.transpose()
            
                gene_pos = gene_variants['Position']
                        
                start = Vcum_array[c]
                end = start + V
                contig_index[start:end] = [gene]*V
                positions[start:end] = gene_pos.as_matrix()        
            except KeyError:
                pass
        return (tauStar,positions,contig_index)
        
    def storeStarState(self,iter):
    
        for c in xrange(self.C):
            if self.gene_ll[c] > self.gene_llstar[c]:
                self.eta_star[c,:] = np.copy(self.eta[c,:])
                self.gene_llstar[c] = self.gene_ll[c]
    