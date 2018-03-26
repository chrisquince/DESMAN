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
import logging

from numpy import array, log, exp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from numpy.random import RandomState

def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring


class Output_Results():

    def __init__(self,outputDir):
        self.outputDir = outputDir
        
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        
        self.log_file_name = self.outputDir+"/log_file.txt"
        
        logging.basicConfig(
            filename=self.log_file_name,
            level=logging.INFO,
            filemode='w', # Overwrites old log file
            format='%(asctime)s:%(levelname)s:%(name)s:%(message)s'
            )
        
        logging.info("Results created in {0}".format(
            os.path.abspath(self.outputDir)))

        print("Up and running. Check {0} for progress".format(
            os.path.abspath(self.log_file_name)), file=sys.stderr)
        
        
    def set_Variants(self,variants):
        self.variants = variants
        self.contig_names = variants.index.tolist()
        self.position = variants['Position']
        
    def set_Variant_Filter(self,variantFilter):
        self.variantFilter = variantFilter
        self.filtered_contig_names = []
        self.filtered_position = []
        for i in self.variantFilter.selected_indices:
            self.filtered_contig_names.append(self.contig_names[i])
            self.filtered_position.append(self.position.iloc[i])
        
    def set_haplo_SNP(self,haplo_SNP,genomes):
        self.haplo_SNP = haplo_SNP
        
        #logLL = haplo_SNP.logLikelihood(haplo_SNP.gamma,haplo_SNP.tau,haplo_SNP.eta)
        
        #AIC = 2.0*haplo_SNP.calcK() - 2.0*logLL
        
        #DIC = haplo_SNP.DIC()
        
        meanDev = haplo_SNP.meanDeviance()
        
        fitFile = self.outputDir+"/fit.txt"
        with open(fitFile, "w") as text_file:
            text_file.write("Fit,%d,%d,%f,%f\n"%(genomes,haplo_SNP.G,haplo_SNP.lp_star, meanDev))

        logging.info("Wrote fit stats") 
        
    def outPredFit(self,haplo_SNP,genomes):
                
        meanDev = haplo_SNP.meanDeviance()
        
        fitFile = self.outputDir+"/fitP.txt"
        with open(fitFile, "w") as text_file:
            text_file.write("Fit,%d,%d,%f,%f\n"%(genomes,haplo_SNP.G,haplo_SNP.lp_star, meanDev))

        logging.info("Wrote pred fit stats") 
        
    def output_Filtered_Tau(self,tau):
        tau_res = np.reshape(tau,(self.haplo_SNP.V,self.haplo_SNP.G*4))
        tau_df = p.DataFrame(tau_res,index=self.filtered_contig_names)
    
        tau_df['Position'] = self.filtered_position
        cols = tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        tau_df = tau_df[cols]
        tau_df.to_csv(self.outputDir+"/Filtered_Tau_star.csv")
        logging.info("Wrote filtered tau star haplotype predictions")
    
    def output_Tau_Mean(self,tauProb):
        
        tau_res = np.reshape(tauProb,(self.haplo_SNP.V,self.haplo_SNP.G*4))
        tau_df = p.DataFrame(tau_res,index=self.filtered_contig_names)
    
        tau_df['Position'] = self.filtered_position
        cols = tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        tau_df = tau_df[cols]
        tau_df.to_csv(self.outputDir+"/Tau_Mean.csv")
        logging.info("Wrote probabilistic tau haplotype predictions")
    
    def output_collated_Tau(self,haplo_SNP_NS,full_variants):
        VS = haplo_SNP_NS.V + self.haplo_SNP.V
        collateTau = np.zeros((VS,self.haplo_SNP.G,4), dtype=np.int)
        collatePTau = np.zeros((VS,self.haplo_SNP.G,4))
        pTau_NS = haplo_SNP_NS.probabilisticTau()
        pTau = self.haplo_SNP.probabilisticTau()
         
        g = 0
        h = 0
        for v in range(VS):
        
            if self.variantFilter.selected[v] != True:
                collateTau[v,:] = haplo_SNP_NS.tau_star[g,:]
                collatePTau[v,:] = pTau_NS[g,:]
                g = g+1
            else:
                collateTau[v,:] = self.haplo_SNP.tau_star[h,:]
                collatePTau[v,:] = pTau[h,:]
                h = h + 1
                
        
        full_contig_names = full_variants.index.tolist()
        full_position = full_variants['Position']
        
        original_contig_names = []
        original_position = []
        for i in self.variantFilter.selected_indices_original:
            original_contig_names.append(full_contig_names[i])
            original_position.append(full_position[i])
        
        collateTau_res = np.reshape(collateTau,(VS,self.haplo_SNP.G*4))
        collatePTau_res = np.reshape(collatePTau,(VS,self.haplo_SNP.G*4))
        
        collate_tau_df = p.DataFrame(collateTau_res,index=original_contig_names)
        collate_tau_df['Position'] = original_position
        
        cols = collate_tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        collate_tau_df = collate_tau_df[cols]
        collate_tau_df.to_csv(self.outputDir+"/Collated_Tau_star.csv")
        logging.info("Wrote all tau haplotype predictions")
        collate_ptau_df = p.DataFrame(collatePTau_res,index=original_contig_names)
        collate_ptau_df['Position'] = original_position
        
        cols = collate_ptau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        collate_ptau_df = collate_ptau_df[cols]
        collate_ptau_df.to_csv(self.outputDir+"/Collated_Tau_mean.csv")
        logging.info("Wrote all probabilistic tau haplotype predictions")    
    
    def output_Gamma_Mean(self,gamma):
        #output max posterior relative frequencies gamma
        
        varCols = self.variants.columns.values.tolist()
        originalS = (len(varCols) - 1)//4
        idx = list(range(1,originalS*4,4))
        sampleNames = [varCols[i] for i in idx] 
        sampleNames = [ rchop(x,'-A') for x in sampleNames ]
        sampleNames = [sampleNames[i] for i in self.variantFilter.sample_indices]
        
        gamma_df = p.DataFrame(gamma,index=sampleNames)
        gamma_df.to_csv(self.outputDir+"/Gamma_mean.csv")
        logging.info("Wrote mean gamma haplotype relative frequencies")
    
    def output_Gamma(self,gamma):
        #output max posterior relative frequencies gamma
        
        varCols = self.variants.columns.values.tolist()
        originalS = (len(varCols) - 1)//4
        idx = list(range(1,originalS*4,4))
        sampleNames = [varCols[i] for i in idx] 
        sampleNames = [ rchop(x,'-A') for x in sampleNames ]
        sampleNames = [sampleNames[i] for i in self.variantFilter.sample_indices]
        
        gamma_df = p.DataFrame(gamma,index=sampleNames)
        gamma_df.to_csv(self.outputDir+"/Gamma_star.csv")
        logging.info("Wrote gamma haplotype relative frequencies")
    
    
    def output_Eta(self,eta):
    
        eta_df = p.DataFrame(eta)
        eta_df.to_csv(self.outputDir+"/Eta_star.csv")
        logging.info("Wrote transition error matrix")
    
    def output_Eta_Mean(self,eta):
    
        eta_df = p.DataFrame(eta)
        eta_df.to_csv(self.outputDir+"/Eta_mean.csv")
        logging.info("Wrote transition error matrix")
    
    
    def output_Selected_Variants(self):
        #write out selected variants really need to remove filtered samples here tooo...
        selected_Variants = self.variants[self.variantFilter.selected]
        selected_Variants.to_csv(self.outputDir+"/Selected_variants.csv")
        logging.info("Wrote selected variants")
    def output_Pickled_haploSNP(self):
        
        with open(self.outputDir+"/haplo_SNP.pickle", 'w') as f:
            pickle.dump(self.haplo_SNP,f)
        logging.info("Wrote pickled haplo_SNP object")