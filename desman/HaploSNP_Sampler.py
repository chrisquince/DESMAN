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

#user defined modules
import Variant_Filter as vf
import Init_NMFT as inmft
import Desman_Utils as du

class Constants(object):
    MAX_LOG_DIR_PROB = 100.0

class HaploSNP_Sampler():
    
    def __init__(self,snps,G,randomState,fixed_tau=None,burn_iter=None,max_iter=None,alpha_constant=10.0,delta_constant=10.0, epsilon=1.0e-6):

        if burn_iter is None:
            self.burn_iter = 250
        else:
            self.burn_iter = burn_iter

        if max_iter is None:
            self.max_iter = 250
        else:
            self.max_iter = max_iter
            
        self.tau_comp_iter = 10

        self.randomState = randomState
        self.G = G

        self.V = snps.shape[0] #number of variants
        self.S = snps.shape[1]

        #set read counts per contig per sample
        self.variants = snps 
        
        self.epsilon = epsilon
        
        #create matrix of genome frequencies gamma and initialise in each site from Dirichlet
        self.delta = np.empty(4); self.delta.fill(delta_constant)
        self.delta_constant = delta_constant
        
        self.alpha = np.empty(self.G); self.alpha.fill(alpha_constant)
        self.alpha_constant = alpha_constant
        
        self.gamma = self.randomState.dirichlet(self.alpha, size=self.S)
        self.gamma_store = np.zeros((self.max_iter,self.S,self.G))
        
        #assignments of genomes to SNPs
        if fixed_tau is None: 
            self.tau = np.zeros((self.V,self.G,4), dtype=np.int)
            self.tau_store = np.zeros((self.max_iter,self.V,self.G,4), dtype=np.int)
            self.tauIndices = np.zeros((self.V),dtype=np.int)
            #assign randomly
            tri = self.randomState.randint(0, 4, self.V*self.G)
            trir = np.reshape(tri,(self.V,self.G))
            for v in range(self.V):
                for g in range(self.G):
                    self.tau[v,g,trir[v,g]] = 1
        else:
            self.tau = np.reshape(fixed_tau,(self.V,self.G,4))
            self.tau_store = np.zeros((self.max_iter,self.V,self.G,4), dtype=np.int)
            self.tauIndices = np.zeros((self.V),dtype=np.int)
            
        
        #initial error transition matrix rate
        self.eta = 0.96*np.identity((4)) + 0.01*np.ones((4,4))
        self.eta_store = np.zeros((self.max_iter,4,4))
        
        #numbers of bases of type a deriving from b with VXSXAXB
        self.E = np.zeros((self.V,self.S,4,4),dtype=np.int)
        self.E_store = np.zeros((self.max_iter,self.V,self.S,4,4),dtype=np.int)
        
        #assignment of bases to genomes
        self.mu = np.zeros((self.V,self.S,4,self.G),dtype=np.int)
        self.mu_store = np.zeros((self.max_iter,self.V,self.S,4,self.G),dtype=np.int)
        
        t1 = np.tile(np.arange(4,dtype=np.int),(self.G,1))
        temparray = du.cartesian(t1)
        self.nTauStates = 4 ** self.G;
        #stores all possible assignments of states to genomes TXGX4
        self.tauStates = np.zeros((self.nTauStates,self.G,4),dtype=np.int)
        
        for t in range(self.nTauStates):
            for g in range(self.G):
                self.tauStates[t,g,temparray[t,g]] = 1
                
            
        #useful to store vectorized matrix of base assignments
        self.amatrix = np.identity(4, dtype=np.int)
        
        self.ll = 0.0
        
        self.tauMap = np.zeros((self.G,4),dtype=np.int)
        
        for g in range(G):
            for a in range(4):
                self.tauMap[g,a] = a*(4**(G - g - 1))
                
    
    def calcK(self):
    
        return self.V*self.G + self.S*(self.G - 1)
            
    def sampleLogProb(self,adLogProbS):
        dP = np.exp(adLogProbS - np.max(adLogProbS))
        dP = dP/np.sum(dP,axis=0)
        return np.flatnonzero(self.randomState.multinomial(1,dP,1))[0]
    
    def baseProbabilityGivenTau(self,tauState,gamma,eta):
        #returns 4 dim vector of each base probability A,C,G,T
        #given gamma and defined assignment tauState
        #loop bases

        return np.einsum('jk,lj,km->lm',tauState,gamma,eta)
        #tensortdot(np.dot(tauState,eta),gamma,axes=([1],[1])) 
    
    def tauDist(self,tau1,tau2):
        dist = 0
        for g in range(self.G):
            id1 = np.where(tau1[g,:] == 1)[0]
            id2 = np.where(tau2[g,:] == 1)[0]
            
            if(id1[0] != id2[0]):
                dist += 1
        
        return dist 
    
    def sampleTau(self):

        nchange = 0
        for v in range(self.V):
            #calculate probability of assignment of each genome to 1 of 4 bases
            for g in range(self.G):
                propTau = np.zeros((4,self.G,4),dtype=np.int)    
                stateLogProb = np.zeros(4)
                
                for a in range(4):
                    propTau[a,:,:] = self.tau[v,:,:]
                    propTau[a,g,:] = np.zeros(4,dtype=np.int)
                    propTau[a,g,a] = 1
                        
                    siteProb  = self.baseProbabilityGivenTau(propTau[a,:,:],self.gamma,self.eta)
                    st1 = np.log(siteProb)*self.variants[v,:,:]
                    stateLogProb[a] = st1.sum()
                    
                s = self.sampleLogProb(stateLogProb)
                    
                self.tau[v,:,:] = propTau[s,:,:] 
                    
                tsample = self.mapTauState(self.tau[v,:,:])
            
                if self.tauIndices[v] != tsample:
                    nchange+=1
                    
                    self.tauIndices[v] = tsample
        return nchange
        
    def mapTauState(self,tauState):
        map = np.einsum('ga,ga',self.tauMap,tauState)
        return map

    def updateTauIndices(self):
        for v in range(self.V):
            tidx = self.mapTauState(self.tau[v,:,:])
            self.tauIndices[v] = tidx
            
    def assignTau(self,assignMatrix):
        """Computes tau matrix for new sets of variants"""
        N = assignMatrix.shape[0]
        assignVariants = np.reshape(assignMatrix, (N,self.S,4)) 
        
        assignTau = np.zeros((N,self.G,4), dtype=np.int)
        conf = np.zeros(N)
        #first compute base probabilities at each site for each state
        siteProb = np.zeros((self.nTauStates,self.S,4))
        
        for t in range(self.nTauStates): 
                siteProb[t,:,:] = self.baseProbabilityGivenTau(self.tauStates[t,:,:],self.gamma_star,self.eta_star)
        stateLogProb = np.zeros(self.nTauStates)
        #loop each variant
        for n in range(N):
            #calculate probability of all possible assignments of genomes to variants
            #at each SNP
            #there are 4^G of these
            for t in range(self.nTauStates):
                st1 = np.log(siteProb[t,:,:])*assignVariants[n,:,:]
                
                stateLogProb[t] = st1.sum()
            tsample = self.sampleLogProb(stateLogProb)
            dP = np.exp(stateLogProb - np.max(stateLogProb))
            dP = dP/np.sum(dP,axis=0)
            conf[n] = np.amax(dP)
            #print "Assign " + str(n) + " " + str(conf[n])
            assignTau[n,:,:] = self.tauStates[tsample,:,:]
        return (assignTau,conf)
        
    def sampleGamma(self):
        #sample gamma from Dirichlet in each sample
        #get frequencies SXG sites X genomes
        sum_mu = self.mu.sum(axis=(0,2))
        #loop samples setting genome frequencies
        for s in range(self.S):
            self.gamma[s,:] = self.randomState.dirichlet(self.alpha + sum_mu[s,:])
        
        self.gamma[self.gamma < self.epsilon] = self.epsilon 
        row_sums = self.gamma.sum(axis=1)
        self.gamma = self.gamma / row_sums[:, np.newaxis]
        
    def sampleEta(self):
        Esum =  self.E.sum(axis=(0,1)) #AXB with A deriving from B
        
        #self.E = np.zeros((self.V,self.S,4,4),dtype=np.int)
        
        for a in range(4):
            self.eta[a,:] = self.randomState.dirichlet(self.delta + Esum[:,a])
        
        
    def sampleMu(self,tauC,gammaC,etaC):
        
        tau_gamma_eta = np.einsum('ijk,lj,km->ilmkj',tauC,gammaC,etaC)
  
        tau_gamma_eta_gsum = tau_gamma_eta.sum(axis=4)
  
        #loop each variant
        for v in range(self.V):
            #loop each site
            for s in range(self.S):
                
                trans_matrix = tau_gamma_eta_gsum[v,s,:,:] 
                tsum = trans_matrix.sum(axis=1)
                trans_matrix = trans_matrix/tsum[:,np.newaxis] 
                
                for a in range(4):
                    #AXB with A observed deriving from B true
                    self.E[v,s,a,:] = self.randomState.multinomial(self.variants[v,s,a], trans_matrix[a,:])
                    
                    self.mu[v,s,a,:] = np.zeros((self.G),dtype=np.int)
                    
                    for b in range(4):
                        if(self.E[v,s,a,b] > 0):
                            tge = tau_gamma_eta[v,s,a,b,:]
                            tge = tge/tge.sum()
                            self.mu[v,s,a,:] += self.randomState.multinomial(self.E[v,s,a,b], tge)
                     
                       
                 
    def burn(self): #perform max_iter Gibbs updates
        iter = 0
        while (iter < self.burn_iter):
            self.sampleMu(self.tau,self.gamma,self.eta)
            self.sampleGamma()
            self.sampleTau()
            self.sampleEta()
            self.ll = self.logLikelihood(self.gamma,self.tau,self.eta)
            print str(iter) + " " + str(self.ll)
            
            iter = iter + 1
    
    def storeStarState(self,iter):
        self.gamma_star = np.copy(self.gamma)
        self.tau_star = np.copy(self.tau)
        self.tauIndices_star = np.copy(self.tauIndices)
        self.eta_star = np.copy(self.eta)
        self.iter_star = iter
        self.ll_star = self.ll
    
    def update(self): #perform max_iter Gibbs updates
        iter = 0
        self.ll = self.logLikelihood(self.gamma,self.tau,self.eta)
        self.storeStarState(iter)
        
        while (iter < self.max_iter):
            self.sampleMu(self.tau,self.gamma,self.eta)
            self.sampleGamma()
            
            nchange = self.sampleTau()
                
            self.sampleEta()
            
            self.ll = self.logLikelihood(self.gamma,self.tau,self.eta)
            if(self.ll > self.ll_star):
                self.storeStarState(iter)
            self.mu_store[iter,]=np.copy(self.mu)
            self.tau_store[iter,]=np.copy(self.tau)
            self.E_store[iter,]=np.copy(self.E)
            self.eta_store[iter,] = np.copy(self.eta)
            self.gamma_store[iter,] = np.copy(self.gamma)    
            
            print str(iter) + "," + str(nchange) + "," + str(self.ll)
            sys.stdout.flush()
            iter = iter + 1
    
    def burnTau(self):
        iter = 0
        self.ll = self.logLikelihood(self.gamma_star,self.tau,self.eta_star)
        
        while (iter < self.burn_iter):
            nchange = self.sampleTau()
            
            self.ll = self.logLikelihood(self.gamma_star,self.tau,self.eta_star)
               
            print str(iter) + "," + str(nchange) + "," + str(self.ll)
            sys.stdout.flush()
            iter = iter + 1
    
    
    def updateTau(self): #perform max_iter Gibbs updates
        
        iter = 0
        self.ll = self.logLikelihood(self.gamma_star,self.tau,self.eta_star)
        
        self.ll_star = self.ll
        self.tau_store[iter,]=np.copy(self.tau)
        
        while (iter < self.max_iter):
                    
            nchange = self.sampleTau()
            
            self.ll = self.logLikelihood(self.gamma_star,self.tau,self.eta_star)
            if (self.ll > self.ll_star):
                self.tau_star = np.copy(self.tau)
                self.ll_star = self.ll
             
            self.tau_store[iter,]=np.copy(self.tau)
               
            print str(iter) + "," + str(nchange) + "," + str(self.ll)
            sys.stdout.flush()
            iter = iter + 1
    
    
    def update_fixed_tau(self): #perform max_iter Gibbs updates
        iter = 0
        self.ll = self.logLikelihood(self.gamma,self.tau,self.eta)
        self.storeStarState(iter)
        
        while (iter < self.max_iter):
            self.sampleMu(self.tau,self.gamma,self.eta)
            self.sampleGamma()
            self.sampleEta()
            
            self.ll = self.logLikelihood(self.gamma,self.tau,self.eta)
            if(self.ll > self.ll_star):
                self.storeStarState(iter)
            self.mu_store[iter,]=np.copy(self.mu)
            self.tau_store[iter,]=np.copy(self.tau)
            self.E_store[iter,]=np.copy(self.E)
            self.eta_store[iter,] = np.copy(self.eta)
            self.gamma_store[iter,] = np.copy(self.gamma)    
            #print str(iter) + " " + str(self.ll)
            
            iter = iter + 1
    
    
    def logLikelihood(self,cGamma,cTau,cEta):
        """Computes data log likelihood given parameter states"""
        logLL = 0.0
        
        probVS = np.einsum('ijk,lj,km->ilm',cTau,cGamma,cEta)
        
        #loop each variant
        for v in range(self.V):
            #loop each site
            for s in range(self.S):                    
                logLL += du.log_multinomial_pdf(self.variants[v,s,:], probVS[v,s,:])
        return logLL
    
    def logTauProb(self,cGamma,cEta):
        
        ret = 0.0
    
        #first compute base probabilities at each site for each state
        siteProb = np.zeros((self.nTauStates,self.S,4))
        
        for t in range(self.nTauStates): 
            siteProb[t,:,:] = self.baseProbabilityGivenTau(self.tauStates[t,:,:],cGamma,cEta)
        
        stateLogProb = np.zeros(self.nTauStates)
        #loop each position
        for v in range(self.V):
            #calculate probability of all possible assignments of genomes to variants
            #at each SNP
            #there are 4^G of these
            for t in range(self.nTauStates):
                st1 = np.log(siteProb[t,:,:])*self.variants[v,:,:]
                
                stateLogProb[t] = st1.sum()
            
            dP = np.exp(stateLogProb - np.max(stateLogProb))
            dP = dP/np.sum(dP,axis=0)
        
            ret += log(dP[self.tauIndices_star[v]])
            
        return ret

    
    def chibMarginalLogLikelihood(self):
        #compute likelihood
        cMLogL = self.logLikelihood(self.gamma_star,self.tau_star,self.eta_star)
        
        #add on priors
        for s in range(self.S):
            cMLogL += du.log_dirichlet_pdf(self.gamma_star[s,:], self.alpha)
        
        for a in range(4):
            cMLogL += du.log_dirichlet_pdf(self.eta_star[a,:],self.delta)
        
        #need tau prior assume uniform over all possible states
        for v in range(self.V):
            cMLogL += log(1.0/float(self.nTauStates))
        
        #estimate posteriors
        tauHat = 0.0
        gammaHat = 0.0
        etaHat = 0.0
        
        #compute first tau term
        storeLogTau = np.zeros(self.tau_comp_iter)
        for i in range(self.tau_comp_iter):
            storeLogTau[i] = self.logTauProb(self.gamma_store[i,:],self.eta_store[i,:])
        
        maxSampleTau = np.max(storeLogTau)
        
        storeLogTau = storeLogTau - maxSampleTau
        storeLogTau = np.exp(storeLogTau)
        tauSum = storeLogTau.sum()
        
        logTauHat =  maxSampleTau + np.log(tauSum) - np.log(self.tau_comp_iter) 
        
        #sample for pi term
        storeLogGamma = np.zeros(self.max_iter)
        for i in range(self.max_iter):
                
            self.sampleMu(self.tau_star,self.gamma,self.eta)
            
            self.sampleGamma()
                
            self.sampleEta()
            
            logTotalP = 0.0;
            
            sum_mu = self.mu.sum(axis=(0,2))
            for s in range(self.S):
                logP = du.log_dirichlet_pdf(self.gamma_star[s,:], self.alpha + sum_mu[s,:])
                
                logTotalP += logP
            
            storeLogGamma[i] = logTotalP
            
        maxSampleGamma = np.max(storeLogGamma)
        
        storeLogGamma = storeLogGamma - maxSampleGamma
        storeLogGamma = np.exp(storeLogGamma)
        gammaSum = storeLogGamma.sum()
        
        logGammaHat =  maxSampleGamma + np.log(gammaSum) - np.log(self.max_iter) 
        
        #sample for epsilon term
        storeLogEpsilon = np.zeros(self.max_iter)
        for i in range(self.max_iter):
                
            self.sampleMu(self.tau_star,self.gamma_star,self.eta)
                
            self.sampleEta()
            
            sum_E =  self.E.sum(axis=(0,1))
            
            logTotalE = 0.0
            for a in range(4):
                logE = du.log_dirichlet_pdf(self.eta_star[a,:],self.delta + sum_E[:,a])
                logTotalE += logE
            
            storeLogEpsilon[i] = logTotalE
        
        maxSampleEpsilon = np.max(storeLogEpsilon)
        
        storeLogEpsilon = storeLogEpsilon - maxSampleEpsilon
        storeLogEpsilon = np.exp(storeLogEpsilon)
        epsilonSum = storeLogEpsilon.sum()
        
        logEpsilonHat =  maxSampleEpsilon + np.log(epsilonSum) - np.log(self.max_iter)
        
        print str(cMLogL) +  "," + str(logGammaHat) + "," + str(logEpsilonHat) + "," + str(logTauHat)
        cMLogL += -logGammaHat - logEpsilonHat - logTauHat
        
        return cMLogL
        
    def calculateSND(self, tau):
        """Calculates number of single nucleotide differences between strains given tau"""
        snd = np.zeros((self.G,self.G),dtype=np.int)
        N = tau.shape[0]
        for g in range(self.G):
            snd[g,g] = 0
            for h in range(g):
                overlap = 0.0;
                for v in range(N):
                    idg = np.where(tau[v,g,:] == 1)[0]
                    idh = np.where(tau[v,h,:] == 1)[0]
                    if(idg[0] == idh[0]):
                        overlap += 1 
                
                snd[g,h] = N - overlap
                
                snd[h,g] = snd[g,h]
        
        return snd
        
    def variableTau(self, tau):
        """Calculates positions with variable bases"""
        N = tau.shape[0]
        variable_tau = np.zeros((N), dtype=bool)
        for v in range(N):
            diff = False
            id0 = np.where(tau[v,0,:] == 1)[0]
            for g in range(1,self.G):
                idg = np.where(tau[v,g,:] == 1)[0]
                if(idg[0] != id0[0]):
                    diff = True 
            variable_tau[v] = diff
            
        return variable_tau
        
    def compSND(self, tau1,tau2):
        """Calculates number of single nucleotide differences between strains given tau"""
        
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
                
                #snd[h,g] = snd[g,h]
        
        return snd
    
    def removeDegenerate(self):
        snd = self.calculateSND(self.tau)
        NU = 0
        deleted = np.zeros(self.G,dtype=bool)
        allmapped = []
        for g in range(self.G):
            if not deleted[g]:
                NU = NU + 1
            gmap = []
            
            for h in range(g+1,self.G):
                if not deleted[h]:
                    if snd[g,h] == 0:
                        deleted[h] = True
                        gmap.append(h)
        
            allmapped.append(gmap)
        NU = self.G - deleted.sum()
        
        tau_new = np.zeros((self.V,NU,4), dtype=np.int)
        gamma_new = np.zeros((self.S,NU))
        NU = 0
        for g in range(self.G):
    
            if not deleted[g]:
                tau_new[:,NU,:] = self.tau[:,g,:]
                gamma_new[:,NU] = self.gamma[:,g]
               
                for h in allmapped[g]:
                    gamma_new[:,NU] += self.gamma[:,h]
                
                NU = NU + 1
                
        self.gamma = gamma_new
        self.tau = tau_new
        self.G = NU
        
        self.alpha = np.empty(self.G); self.alpha.fill(self.alpha_constant)
        self.gamma_store = np.zeros((self.max_iter,self.S,self.G))
        self.tau_store = np.zeros((self.max_iter,self.V,self.G,4), dtype=np.int)
        
        #assignment of bases to genomes
        self.mu = np.zeros((self.V,self.S,4,self.G),dtype=np.int)
        self.mu_store = np.zeros((self.max_iter,self.V,self.S,4,self.G),dtype=np.int)
        
        t1 = np.tile(np.arange(4,dtype=np.int),(self.G,1))
        temparray = du.cartesian(t1)
        self.nTauStates = 4 ** self.G;
        #stores all possible assignments of states to genomes TXGX4
        self.tauStates = np.zeros((self.nTauStates,self.G,4),dtype=np.int)
        
        for t in range(self.nTauStates):
            for g in range(self.G):
                self.tauStates[t,g,temparray[t,g]] = 1    
        
        self.tauMap = np.zeros((self.G,4),dtype=np.int)
        
        for g in range(self.G):
            for a in range(4):
                self.tauMap[g,a] = a*(4**(self.G - g - 1))    
    
        self.updateTauIndices()
    
    def probabilisticTau(self):    
        
        sumTau = self.tau_store.sum(axis = 0)
        
        probTau = sumTau/float(self.max_iter)
        
        return probTau