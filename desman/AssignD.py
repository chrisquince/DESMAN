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

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("gamma_file", help="gamma file")

    parser.add_argument("eta_file", help="eta file")

    parser.add_argument("assign_file", help="assign file")

    parser.add_argument('-o','--output_dir', type=str, default="output",
        help=("string specifying output directory and file stubs"))

    #get command line arguments  
    args = parser.parse_args()
    
    import ipdb; ipdb.set_trace()
    
    gamma_file = args.gamma_file
    eta_file = args.eta_file
    assign_file = args.assign_file
    output_dir = args.output_dir
    
    variants    = p.read_csv(assign_file, header=0, index_col=0)
    variants_matrix = variants.as_matrix()
    variants_matrix = np.delete(variants_matrix, 0, 1)
    
    snps = np.reshape(variants_matrix, (variants_matrix.shape[0],variants_matrix.shape[1]/4,4))
    
    #create new random state with fixed seed 
    prng = RandomState(238329)

    gamma_df = p.read_csv(gamma_file, header=0, index_col=0) 
    
    gmatrix = gamma_df.as_matrix()
    G = gmatrix.shape[1]
    
    haplo_SNP = hsnp.HaploSNP_Sampler(snps,G,prng)
    
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
     #  import ipdb; ipdb.set_trace()
        (assignTau,confTau) = haplo_SNP.assignTau(assigns_matrix)
        snda = haplo_SNP.calculateSND(assignTau)
        
        assign_contig_names = assigns.index.tolist()
        assign_position = assigns['Position']
        
        AV = assigns_matrix.shape[0]
        assign_tau_res = np.reshape(assignTau,(AV,haplo_SNP.G*4))
        assign_tau_df = p.DataFrame(assign_tau_res,index=assign_contig_names)
        conf_tau_df = p.DataFrame(confTau,index=assign_contig_names)
        
        assign_tau_df['Position'] = assign_position
        conf_tau_df['Position'] = assign_position
        
        cols = assign_tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        
        assign_tau_df = assign_tau_df[cols]
        assign_tau_df.to_csv(output_dir+"/Assigned_Tau_star.csv")
    
        cols = conf_tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        
        conf_tau_df = conf_tau_df[cols]
        conf_tau_df.to_csv(output_dir+"/Assigned_Tau_conf.csv")
    
if __name__ == "__main__":
    main(sys.argv[1:])
