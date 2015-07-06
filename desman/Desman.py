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
    parser.add_argument("variant_file", help="input SNP frequencies")

    parser.add_argument('-g','--genomes', type=int, required=True,
        help=("specify the genome number"))
        
    parser.add_argument('-f','--filter_variants',nargs='?', const=3.84, type=float, 
        help=("filters variants by negative binomial loge likelihood defaults to 40.0"))

    parser.add_argument('-r','--random_select',nargs='?', const=1e3, type=int, 
        help=("selects subset of variants passing filter to build model and assigns others"))

    parser.add_argument('-e','--eta_file', type=file, 
        help=("reads initial eta matrix from file"))

    parser.add_argument('-a','--assign_file', type=file, 
        help=("calculates haplotype profiles for these SNPs using fitted gamma, eta values"))
    
    parser.add_argument('-o','--output_dir', type=str, default="output",
        help=("string specifying output directory and file stubs"))
    
    parser.add_argument('-p', '--optimiseP', default=True, type=bool,
        help=("optimise proportions in likelihood ratio test"))
        
    parser.add_argument('-l', '--xtraslow', action='store_true',
        help=("always perform full state sampling"))
    
    parser.add_argument('-m','--min_coverage', type=str, default=5.0,
        help=("minimum coverage for sample to be included"))
    
    parser.add_argument('-t','--tau_file', type=str, default=None,
        help=("read in known tau file"))

    parser.add_argument('-q','--max_qvalue',default=1.0e-3, type=float, 
        help=("specifies q value cut-off for variant detection defaults 1.0e-3"))
    
    parser.add_argument('-s','--random_seed',default=23724839, type=int, 
        help=("specifies seed for numpy random number generator defaults to 23724839 applied after random filtering"))
    
    #get command line arguments  
    args = parser.parse_args()
    variant_file = args.variant_file
    assign_file = args.assign_file
    eta_file = args.eta_file
    filter_variants = args.filter_variants
    random_select = args.random_select
    output_dir = args.output_dir
    tau_file = args.tau_file
    min_coverage = args.min_coverage
    optimiseP = args.optimiseP
    max_qvalue = args.max_qvalue
    genomes = args.genomes
    xtraslow = args.xtraslow
    random_seed = args.random_seed
    
    #create new random state with fixed seed 
    prng = RandomState(238329)
    
    #read in snp variants
    variants    = p.read_csv(variant_file, header=0, index_col=0)
    
    import ipdb; ipdb.set_trace()
    
    #construct variant filter to only select most likely SNPS
    variant_Filter = vf.Variant_Filter(variants, randomState = prng, optimise = optimiseP, threshold = filter_variants, min_coverage = min_coverage, qvalue_cutoff = max_qvalue)
    
    if filter_variants is not None:    
        #variants_matrix = variant_Filter.get_filtered_Variants()
        variant_Filter.get_filtered_VariantsLogRatio()
    
    if eta_file is not None:
        eta_df = p.read_csv(eta_file, header=0, index_col=0)    
        variant_Filter.eta = eta_df.as_matrix()
         
    if random_select is not None:
        variant_Filter.select_Random(random_select)
    
    if tau_file is not None:
        tau = p.read_csv(tau_file, header=0, index_col=0)
        tau_matrix = tau.as_matrix()
        tau_matrix = np.delete(tau_matrix, 0, 1)
        haplo_SNP_fixed = HaploSNP_Sampler(snps_filter,tau_matrix.shape[1]/4,tau_matrix)
        haplo_SNP_fixed.update_fixed_tau()
    
    prng = RandomState(args.random_seed)
    init_NMFT = inmft.Init_NMFT(variant_Filter.snps_filter,genomes,prng)
    init_NMFT.factorize()
    
    haplo_SNP = hsnp.HaploSNP_Sampler(variant_Filter.snps_filter,genomes,prng,bSlow=xtraslow)
    haplo_SNP.tau = init_NMFT.get_tau()
    haplo_SNP.updateTauIndices()
    haplo_SNP.gamma = init_NMFT.get_gamma()
    haplo_SNP.eta = variant_Filter.eta
     
    #haplo_SNP.burn()
    haplo_SNP.update()
    #after burn-in phase remove degeneracies?
    haplo_SNP.removeDegenerate()
    
    haplo_SNP.update()
    snd = haplo_SNP.calculateSND(haplo_SNP.tau_star)
    
    if tau_file is not None:
        snd_comp = haplo_SNP.compSND(haplo_SNP_fixed.tau,haplo_SNP.tau)
    #import ipdb; ipdb.set_trace()
    cMLogL = haplo_SNP.chibbMarginalLogLikelihood()
    logLL = haplo_SNP.logLikelihood(haplo_SNP.gamma,haplo_SNP.tau,haplo_SNP.eta)
    snd = haplo_SNP.calculateSND(haplo_SNP.tau_star)
    variableTau = haplo_SNP.variableTau(haplo_SNP.tau_star)
    
    probTau = haplo_SNP.probabilisticTau()
    
    #output log likelihood and Chib's approxn to the log marginal likelihood
    print str(genomes) + "," + str(haplo_SNP.G) + "," + str(logLL) + "," + str(cMLogL)
    
    #output results to files
    #import ipdb; ipdb.set_trace()
    output_Results = outr.Output_Results(variants,haplo_SNP,variant_Filter, output_dir)
    output_Results.output_Filtered_Tau(haplo_SNP.tau_star)
    output_Results.output_Gamma(haplo_SNP.gamma_star)
    output_Results.output_Eta(haplo_SNP.eta_star)
    output_Results.output_Selected_Variants()
    output_Results.output_Pickled_haploSNP()
    
    #assign if assignment file given
    if(assign_file != None):
        assigns    = p.read_csv(assign_file, header=0, index_col=0)    
        assigns_matrix = assigns.as_matrix()
        assigns_matrix = np.delete(assigns_matrix, 0, 1)
        #import ipdb; ipdb.set_trace()
        (assignTau,confTau) = haplo_SNP.assignTau(assigns_matrix)
        snda = haplo_SNP.calculateSND(assignTau)
        
        assign_contig_names = assigns.index.tolist()
        assign_position = assigns['Position']
        
        AV = assigns_matrix.shape[0]
        assign_tau_res = np.reshape(assignTau,(AV,haplo_SNP.G*4))
        assign_tau_df = p.DataFrame(assign_tau_res,index=assign_contig_names)
    
        assign_tau_df['Position'] = assign_position
        cols = assign_tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        assign_tau_df = assign_tau_df[cols]
        assign_tau_df.to_csv(output_dir+"/Assigned_Tau_star.csv")
    
if __name__ == "__main__":
    main(sys.argv[1:])
