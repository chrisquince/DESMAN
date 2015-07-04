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


class Output_Results():

    def __init__(self,variants,haplo_SNP,variantFilter,outputDir):
        self.variantFilter = variantFilter
        self.outputDir = outputDir
        self.variants = variants
        
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        
        self.contig_names = variants.index.tolist()
        self.position = variants['Position']
        self.haplo_SNP = haplo_SNP
        
    def output_Filtered_Tau(self,tau):

        filtered_contig_names = []
        filtered_position = []
        for i in self.variantFilter.selected_indices:
            filtered_contig_names.append(self.contig_names[i])
            filtered_position.append(self.position[i]) 
        
        tau_res = np.reshape(tau,(self.haplo_SNP.V,self.haplo_SNP.G*4))
        tau_df = p.DataFrame(tau_res,index=filtered_contig_names)
    
        tau_df['Position'] = filtered_position
        cols = tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        tau_df = tau_df[cols]
        tau_df.to_csv(self.outputDir+"/Filtered_Tau_star.csv")
    
    def output_Gamma(self,gamma):
        #output max posterior relative frequencies gamma
        
        varCols = self.variants.columns.values.tolist()
        originalS = (len(varCols) - 1)/4
        idx = range(1,originalS*4,4)
        sampleNames = [varCols[i] for i in idx] 
        sampleNames = [sampleNames[i] for i in self.variantFilter.sample_indices]
        
        gamma_df = p.DataFrame(gamma,index=sampleNames)
        gamma_df.to_csv(self.outputDir+"/Gamma_star.csv")
        
    def output_Eta(self,eta):
    
        eta_df = p.DataFrame(eta)
        eta_df.to_csv(self.outputDir+"/Eta_star.csv")
    
    def output_Selected_Variants(self):
        #write out selected variants really need to remove filtered samples here tooo...
        selected_Variants = self.variants[self.variantFilter.selected]
        selected_Variants.to_csv(self.outputDir+"/Selected_variants.csv")
        
    def output_Pickled_haploSNP(self):
        
        with open(self.outputDir+"/haplo_SNP.pickle", 'w') as f:
            cPickle.dump(self.haplo_SNP,f)