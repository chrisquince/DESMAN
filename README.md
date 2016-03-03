# DESMAN _De novo_ Extraction of Strains from MetAgeNomes

![alt tag](desmans.jpg)

##Installation

To install simply type:
    
    sudo python ./setup.py install
    
These items are prerequisities for the installation of desman:

*python v2.7.*
*gcc
*gsl

The installation procedure varies on different systems, 
and described in this README is only how to proceed with a linux (ubuntu) distribution.

The first item, python v2.7.*, should be installed on a modern Ubuntu distribution. 
A c-compiler, e.g. gcc, is needed to compile the c parts of concoct that uses the 
GNU Scientific Library gsl. For linux (ubuntu) this is installed through:

    sudo apt-get install build-essential libgsl0-dev

##Simple example

To illustrate the actual strain inference algorithm we will start with a sample example using base frequencies 
that have been prepared. The starting point for a Desman analysis is a csv file with base frequencies e.g. 

[Strain mock community frequencies for COG0015](data/contig_6or16_genesL_scgCOG0015.freq)

This has a simple format, a header line followed by one line for each position with format:

gene name, position, freq. of A in sample 1, freq. of C in 1,freq. of G in 1,freq. of T in 1,..., freq. of A in sample N, freq. of C in N,freq. of G in N,freq. of T in N 

##Complete _de novo_ strain level analysis from metagenome data


First step is to identify variant positions. This is performed by the program Variant_Filter.py. 
Start assuming you are in the DESMAN repo directory by making a test folder.

    mkdir test
    cd test

Then run the example data file which corresponds to a single COG from the mock community data set 
described in the manuscript. This COG0015 has 933 variant positions. The input file is in the data 
folder. We run the variant filtering as follows:

    python ../desman/Variant_Filter.py data/contig_6or16_genesL_scgCOG0015.freq -o COG0015_out -p

The '-p' flag uses one dimenisonal optimisition to find individual base frequencies.

This version now uses Gibbs sampling for tau.

This generates output files: 

