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
import copy

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
import sampletau

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
    
    parser.add_argument('-i','--no_iter',nargs='?', const=250, type=int, 
        help=("Number of iterations of Gibbs sampler"))
    
    parser.add_argument('-m','--min_coverage', type=float, default=5.0,
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
    no_iter = args.no_iter
    random_seed = args.random_seed
    
    #create new random state with fixed seed 
    prng = RandomState(238329)
    
    #read in snp variants
    variants    = p.read_csv(variant_file, header=0, index_col=0)
    
    #import ipdb; ipdb.set_trace()
    
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
    sampletau.initRNG()
    sampletau.setRNG(args.random_seed)
    
    init_NMFT = inmft.Init_NMFT(variant_Filter.snps_filter,genomes,prng)
    init_NMFT.factorize()
    
    haplo_SNP = hsnp.HaploSNP_Sampler(variant_Filter.snps_filter,genomes,prng,max_iter=no_iter)
    
    haplo_SNP.tau = np.copy(init_NMFT.get_tau(),order='C')
    haplo_SNP.updateTauIndices()
    haplo_SNP.gamma = np.copy(init_NMFT.get_gamma(),order='C')
    haplo_SNP.eta = np.copy(variant_Filter.eta,order='C')
     
    #haplo_SNP.burn()
    
    haplo_SNP.update()
    #after burn-in phase remove degeneracies?
    haplo_SNP.removeDegenerate()
    
    haplo_SNP.update()
    snd = haplo_SNP.calculateSND(haplo_SNP.tau_star)
    
    if tau_file is not None:
        snd_comp = haplo_SNP.compSND(haplo_SNP_fixed.tau,haplo_SNP.tau)
    #import ipdb; ipdb.set_trace()
    #cMLogL = haplo_SNP.chibMarginalLogLikelihood()
    logLL = haplo_SNP.logLikelihood(haplo_SNP.gamma,haplo_SNP.tau,haplo_SNP.eta)
    snd = haplo_SNP.calculateSND(haplo_SNP.tau_star)
    variableTau = haplo_SNP.variableTau(haplo_SNP.tau_star)
    
    probTau = haplo_SNP.probabilisticTau()
    
    #output log likelihood and Chib's approxn to the log marginal likelihood
    #import ipdb; ipdb.set_trace()
    #chib = haplo_SNP.chibMarginalLogLikelihood2()
    
    AIC = 2.0*haplo_SNP.calcK() - 2.0*logLL
    print "Fit" + "," + str(genomes) + "," + str(haplo_SNP.G) + "," + str(logLL) + "," + str(AIC) 
    
    #output results to files
    output_Results = outr.Output_Results(variants,haplo_SNP,variant_Filter, output_dir)
    
    output_Results.output_Filtered_Tau(haplo_SNP.tau_star)

    output_Results.output_Prob_Tau(probTau)
    
    output_Results.output_Gamma(haplo_SNP.gamma_star)
    
    output_Results.output_Eta(haplo_SNP.eta_star)
    
    output_Results.output_Selected_Variants()
    
    #import ipdb; ipdb.set_trace()
    if random_select is not None:
        VS = variant_Filter.snps_filter_original.shape[0] 
        
        snps_notselected = variant_Filter.snps_filter_original[variant_Filter.selected != True,:]
    
        init_NMFT_NS = inmft.Init_NMFT(snps_notselected,haplo_SNP.G,haplo_SNP.randomState)
        
        init_NMFT_NS.gamma = np.transpose(haplo_SNP.gamma)
        
        init_NMFT_NS.factorize_tau()
        
        haplo_SNP_NS = hsnp.HaploSNP_Sampler(snps_notselected,haplo_SNP.G,haplo_SNP.randomState,max_iter=no_iter)
    
        haplo_SNP_NS.tau = init_NMFT_NS.get_tau()
        haplo_SNP_NS.updateTauIndices()
        haplo_SNP_NS.gamma_star = np.copy(haplo_SNP.gamma_star,order='C')
        haplo_SNP_NS.eta_star = np.copy(haplo_SNP.eta_star,order='C')     
        
        haplo_SNP_NS.updateTau()
    
        haplo_SNP_NS.updateTau()

        collateTau = np.zeros((VS,haplo_SNP.G,4), dtype=np.int)
        collatePTau = np.zeros((VS,haplo_SNP.G,4))
        pTau_NS = haplo_SNP_NS.probabilisticTau()
        pTau = haplo_SNP_NS.probabilisticTau()
         
        g = 0
        h = 0
        for v in range(VS):
        
            if variant_Filter.selected[v] != True:
                collateTau[v,:] = haplo_SNP_NS.tau_star[g,:]
                collatePTau[v,:] = pTau_NS[g,:]
                g = g+1
            else:
                collateTau[v,:] = haplo_SNP.tau_star[h,:]
                collatePTau[v,:] = pTau[h,:]
                h = h + 1
                
        
        contig_names = variants.index.tolist()
        position = variants['Position']
        
        original_contig_names = []
        original_position = []
        for i in variant_Filter.selected_indices_original:
            original_contig_names.append(contig_names[i])
            original_position.append(position[i])
        
        collateTau_res = np.reshape(collateTau,(VS,haplo_SNP.G*4))
        collatePTau_res = np.reshape(collatePTau,(VS,haplo_SNP.G*4))
        
        collate_tau_df = p.DataFrame(collateTau_res,index=original_contig_names)
        collate_tau_df['Position'] = original_position
        
        cols = collate_tau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        collate_tau_df = collate_tau_df[cols]
        collate_tau_df.to_csv(output_dir+"/Collated_Tau_star.csv")
        
        collate_ptau_df = p.DataFrame(collatePTau_res,index=original_contig_names)
        collate_ptau_df['Position'] = original_position
        
        cols = collate_ptau_df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        collate_ptau_df = collate_ptau_df[cols]
        collate_ptau_df.to_csv(output_dir+"/Collated_PTau_star.csv")
    #assign if assignment file given
    if(assign_file != None):
        assigns    = p.read_csv(assign_file, header=0, index_col=0)    
        assigns_matrix = assigns.as_matrix()
        assigns_matrix = np.delete(assigns_matrix, 0, 1)
        import ipdb; ipdb.set_trace()
       
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
    
    sampletau.freeRNG()
if __name__ == "__main__":
    main(sys.argv[1:])
