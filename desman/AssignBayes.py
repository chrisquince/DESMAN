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
import HaploSNP_Sampler as hsnp
import Output_Results as outr

def writeTauFile(tauMatrix,fileName,positions,contig_names,G):
    
    AV = tauMatrix.shape[0]
    tauMatrix_res = np.reshape(tauMatrix,(AV,G*4))
    
    tauMatrix_df = p.DataFrame(tauMatrix_res,index=contig_names)
    
    tauMatrix_df['Position'] = positions
    cols = tauMatrix_df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
        
    tauMatrix_df = tauMatrix_df[cols]
    tauMatrix_df.to_csv(fileName)

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gamma_file", help="gamma file")

    parser.add_argument("eta_file", help="eta file")

    parser.add_argument("assign_file", help="assign file")

    parser.add_argument('--fast', dest='fast', action='store_true')

    parser.add_argument('-o','--output_dir', type=str, default="output",
        help=("string specifying output directory and file stubs"))

    parser.add_argument('-i','--max_iter', type=int, default=100,
        help=("string specifying output directory and file stubs"))

    #get command line arguments  
    args = parser.parse_args()
    
    #import ipdb; ipdb.set_trace()
    
    gamma_file = args.gamma_file
    eta_file = args.eta_file
    assign_file = args.assign_file
    output_dir = args.output_dir
    fast = args.fast
    max_iter = args.max_iter
    
    variants    = p.read_csv(assign_file, header=0, index_col=0)
    variants_matrix = variants.as_matrix()
    variants_matrix = np.delete(variants_matrix, 0, 1)
    
    snps = np.reshape(variants_matrix, (variants_matrix.shape[0],variants_matrix.shape[1]/4,4))
    
    #create new random state with fixed seed 
    prng = RandomState(238329)

    gamma_df = p.read_csv(gamma_file, header=0, index_col=0) 
    
    gmatrix = gamma_df.as_matrix()
    G = gmatrix.shape[1]
    
    haplo_SNP = hsnp.HaploSNP_Sampler(snps,G,prng,burn_iter=20,max_iter=max_iter)
    
    #haplo_SNP.tau_star = haplo_SNP.tau 
    
    eta_df = p.read_csv(eta_file, header=0, index_col=0)    
    haplo_SNP.eta_star = eta_df.as_matrix()
    
   
    haplo_SNP.gamma_star = gamma_df.as_matrix()
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    if(assign_file != None):
        assigns    = p.read_csv(assign_file, header=0, index_col=0)    
        assigns_matrix = assigns.as_matrix()
        assigns_matrix = np.delete(assigns_matrix, 0, 1)
        
        haplo_SNP.burnTau()
        
        haplo_SNP.updateTau()
        
        snda = haplo_SNP.calculateSND(haplo_SNP.tau_star)
        
        fileNameS=output_dir+"/Assigned_Tau_star.csv"
        
        assign_positions = assigns['Position']
        assign_contig_names = assigns.index.tolist()
        
        writeTauFile(haplo_SNP.tau_star,fileNameS,assign_positions,assign_contig_names,haplo_SNP.G)
        
        probTau = haplo_SNP.probabilisticTau()
        
        fileNameP=output_dir+"/Assigned_Tau_prob.csv"
        writeTauFile(probTau, fileNameP, assign_positions, assign_contig_names,haplo_SNP.G)
    
        for g in range(haplo_SNP.G):
            uconf = 1.0 - np.max(probTau[:,g,:],axis=1)
            
            avuconf = uconf.sum()/int(haplo_SNP.V)
            
            print str(g)+","+str(avuconf)
if __name__ == "__main__":
    main(sys.argv[1:])
