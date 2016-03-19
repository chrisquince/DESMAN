# DESMAN _De novo_ Extraction of Strains from MetAgeNomes

![alt tag](desmans.jpg)

##Installation

To install simply type:
    
    sudo python ./setup.py install
    
These items are prerequisities for the installation of desman:

1. *python v2.7.*
2. *gcc
3. *gsl

The installation procedure varies on different systems, 
and described in this README is only how to proceed with a linux (ubuntu) distribution.

The first item, python v2.7.*, should be installed on a modern Ubuntu distribution. 
A c-compiler, e.g. gcc, is needed to compile the c parts of concoct that uses the 
GNU Scientific Library gsl. For linux (ubuntu) this is installed through:

```
    sudo apt-get install build-essential libgsl0-dev
```

For convenience we also recommend adding the scripts directory to your path:

```
export PATH=$HOME/myinstalldir/DESMAN/scripts:$PATH

````

Obviously replacing myinstalldir as appropriate and adding this to your .bash_profile file.

##Simple example

To illustrate the actual strain inference algorithm we will start with a simple example using base frequencies 
that have been pre-prepared. Below we also give [a complete example](#complete_example) including 
pre-processing. The starting point for a Desman analysis is a csv file with base frequencies e.g.: 

[Strain mock community frequencies for COG0015](data/contig_6or16_genesL_scgCOG0015.freq)

This has the following format:

    Contig,Position,SampleName1-A,SampleName1-C,SampleName1-G,SampleName1-T,...,SampleNameN-A,SampleNameN-C,SampleNameN-G,SampleNameN-T

where SampleName1,...,SampleNameN gives the names of the different samples in the analysis. Followed 
by one line for each position with format:

    gene name, position, freq. of A in sample 1, freq. of C in 1,freq. of G in 1,freq. of T in 1,..., freq. of A in sample N, freq. of C in N,freq. of G in N,freq. of T in N 


###Finding variant positions for the test data set

The first step is to identify variant positions. This is performed by the desman script Variant_Filter.py. 
Start assuming you are in the DESMAN repo directory by making a test folder.

    mkdir test
    cd test

Then run the example data file which corresponds to a single COG from the mock community data set 
described in the manuscript. This COG0015 has 933 variant positions. The input file is in the data 
folder. We run the variant filtering as follows:

    python ../desman/Variant_Filter.py ../data/contig_6or16_genesL_scgCOG0015.freq -o COG0015_out -p

The variant filtering has a number of optional parameters to see them run:

    python ../desman/Variant_Filter.py -h
    
They should all be fairly self explanatory. We recommend always using the 
the '-p' flag for one dimenisonal optimisition of individual base frequencies if it is not 
too time consuming. The '-o' option is a file stub all output files will be generated with this prefix.
A log file will be generated 'COG0015_out_log.txt' and output files: 

1. COG0015_outp_df.csv: This gives p-values for each position.

2. COG0015_outq_df.csv: This gives q-values for each position.

3. COG0015_outr_df.csv: This gives log-ratio statistics for each position.

4. COG0015_outsel_var.csv: This is the file of selected variants.

5. COG0015_outtran_df.csv: A matrix of estimated error rates.

###Inferring haplotypes and abundances for the test data set

Having found the variant positions we will now the run the program for inferring haplotypes and their abundance:

    desman COG0015_outsel_var.csv -g 5 -e COG0015_outtran_df.csv -o COG0015_out_g5 -i 50 

These parameters specify the variants file. Then number of haplotypes as five '-g 5', an initial 
estimate for the error transition matrix taken from the variant detection '-e COG0015_outtran_df.csv', 
an output directory '-o COG0015_out_g5' and the number of iterations, '-i 50'.
The program takes the selected variants and infers haplotypes and their abundances using the Gibbs sampler given 
the assumption that five strains are present. All output files will be generated in the directory COG0015_out_g5. 
Once the program has finished running a few minutes on a typical computer it will generate the following 
files inside the output directory:

1. log_file.txt: This logs the progress of the program through the three stages: 
NTF initialisation, 'burn-in' Gibbs sampler and the sampling itself. 

2. Eta_star.csv: Prediction for error transition matrix (rows true bases, columns observed probabilities) 						

3. Filtered_Tau_star.csv: Prediction for strain haplotypes. Each row of comma separated file contains:

```
    gene name, position, haplotype1-A,  haplotype1-C,  haplotype1-G,  haplotype1-T,..., haplotypeG-A,  haplotypeG-C,  haplotypeG-G,  haplotypeG-T  
```

where 1 indicates the base present in that haplotype at that position.

4. Gamma_star.csv: This gives the relative frequency of each haplotye in each sample. One row for each sample.

5. Probabilistic_Tau.csv: Prediction for strain haplotypes but now an estimate of the poseterior probability of asignments.

6. Selected_variants.csv: Variants used for strain calling if filtering applied.

7. fit.txt: Statistics evaluating fit as number of haplotypes, number of non-degenerate haplotypes inferred, log-Likelihood, Aikake information criterion

#Complete example of _de novo_ strain level analysis from metagenome data
<a name="complete_example"></a>

To provide an in depth illustration of how to use Deman we will give a complete worked example from a subset of the synthetic community 
used in Quince et al. 2016. We have provided 16 samples, subsampled to 1 million reads from the 64 samples with 11.75 million reads used 
originally. This example is therefore more tractable but the following analysis assumes you have access to a multi-core server. 
We also assume that you have some standard and not so standard sequence analysis software installed:

1. [megahit](https://github.com/voutcn/megahit): A highly efficient metagenomics assembler currently our default for most studies

2. [bwa](https://github.com/lh3/bwa): Necessary for mapping reads onto contigs

3. [bam-readcount](https://github.com/genome/bam-readcount): Used to get per sample base frequencies at each position

4. [samtools] (http://www.htslib.org/download/): Utilities for processing mapped files

5. [CONCOCT](https://github.com/BinPro/CONCOCT): Our own contig binning algorithm

Click the link associated with each application for installation details. To begin obtain the reads from Dropbox:

Untar and unzip the example directory and move into it:

```bash
tar -xvzf Example.tar.gz
cd Example
```

Then assemble the reads. We recommend megahit for this:
```bash
nohup megahit -1 $(<R1.csv) -2 $(<R2.csv) -t 36 -o Assembly --presets meta > megahit.out&
```

We will now perform CONCOCT binning of these contigs. As explained in [Alneberg et al.](http://www.nature.com/nmeth/journal/v11/n11/full/nmeth.3103.html) 
there are good reasons to cut up contigs prior to binning. We will use a script from CONCOCT to do this. So make sure the 
CONCOCT scripts directory is in your path:

```bash
export CONCOCT=$HOME/mypathtoConcoct/CONCOCT
```

Then cut up contigs and place in new dir:

```bash
mkdir contigs
python $CONCOCT/scripts/cut_up_fasta.py -c 10000 -o 0 -m Assembly/final.contigs.fa > contigs/final_contigs_c10K.fa
```

Having cut-up the contigs the next step is to map all the reads from each sample back onto them. First index the contigs with bwa:

```bash
cd contigs
bwa index final_contigs_c10K.fa
cd ..
```

Then perform actual mapping you may want to put this in a shell script:

```bash
mkdir Map

for file in *R1.fastq
do 
   
   stub=${file%_R1.fastq}

   echo $stub

   file2=${stub}_R2.fastq

   bwa mem -t 32 contigs/final_contigs_c10K.fa $file $file2 > Map/${stub}.sam
done
```

```bash
python Lengths.py -i contigs/final_contigs_c10K.fa > contigs/final_contigs_c10K.len
```


```bash
for file in Map/*.sam
do
    stub=${file%.sam}
    stub2=${stub#Map\/}
    echo $stub	
    (samtools view -h -b -S $file > ${stub}.bam; samtools view -b -F 4 ${stub}.bam > ${stub}.mapped.bam; samtools sort -m 1000000000 ${stub}.mapped.bam ${stub}.mapped.sorted; bedtools genomecov -ibam ${stub}.mapped.sorted.bam -g contigs/final_contigs_c10K.len > ${stub}_cov.txt)&
done

```
