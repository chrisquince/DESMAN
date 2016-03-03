# DESMAN _De novo_ Extraction of Strains from MetAgeNomes

![alt tag](https://github.com/chrisquince/DESMAN/desmans.jpg)

##Installation

To install simply type:
    python ./set
    

##Simple example

##Complete _de novo_ strain level analysis from metagenome data


compile Cython module run:
    python setup.py build_ext --inplace

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

