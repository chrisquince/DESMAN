# DESMAN complete example using Snakemake

![alt tag](desmans.jpg)

## Table of Contents
[Installation](#installation)

[Taxonomic profiling](#simple_example)

[Complete example](#complete_example)

<a name="installation"/>

##Installation

### Install Dependencies
DESMAN requires some parts of it to be run with python 3 and some with python 2.7, therefore we strongly recommend the use of Conda. Conda is a tool to isolate your python installation, which allows you to have multiple parallel installations using different versions for different packages, and gives you a very convenient and fast way to install the most common scientific python packages. Additionaly you can use Conda to install all kinds of software into an somewhat isolated environment. Conda is free of charge and you can read more about it [here](https://github.com/conda/conda).

To install the recommended minimal version of conda, miniconda:

```
    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh
    chmod +x Miniconda-latest-Linux-x86_64.sh
    ./Miniconda-latest-Linux-x86_64.sh
    export PATH=~/miniconda2/bin:$PATH
```
Select default answers to all questions and yes at each step note that the miniconda2 bin will then 
be permanently added to the .bashrc so the export will not be needed on future sessions.


After installing Conda, create an new environment that will contain the python 2.7 installation called 'desman_python2_env':

```
    conda create -n desman_python2_env python=2.7
```

and one environment that will contain the python 3 installation called 'desman_python3_env':


```
    conda create -n desman_python3_env python=3
```

In the desman_python2_env environment we'll install the CONCOCT and DESMAN dependencies:

```
    source activate desman_python2_env
    conda install cython numpy scipy biopython pandas pip scikit-learn pysam
```

These items are prerequisities for the installation of DESMAN:

1. *gcc
2. *gsl

The installation procedure varies on different systems,
and described in this README is only how to proceed with a linux (ubuntu) distribution.

A c-compiler, e.g. gcc, is needed to compile the c parts of CONCOCT and DESMAN that use the
GNU Scientific Library GSL. For linux (ubuntu) this is installed through:

```
    sudo apt-get install build-essential libgsl0-dev
```

We then install CONCOCT in a directory ~/repos but this could be any where on your system just change 
as appropriate:

```
    mkdir ~/repos
    cd ~/repos
    git clone https://github.com/BinPro/CONCOCT.git
    cd CONCOCT
    python ./setup.py install
```

and the DESMAN software:

```
    cd ~/repos
    git clone https://github.com/chrisquince/DESMAN.git
    python ./setup.py install
```

In the desman_python3_env we'll install a few dependencies that are not available in the default Conda distribution, instead we install these using the bioconda channel.

```
    source activate desman_python3_env
    conda install -c bioconda snakemake bwa samtools prodigal diamond r r-gplots r-getopt
```

Need to look at the above dependencies failing on my fresh Ubuntu. Instead install through repos e.g. BWA:
We also assume that you have some standard and not so standard sequence analysis software installed:

1. [megahit](https://github.com/voutcn/megahit): A highly efficient metagenomics assembler currently our default for most studies

2. [bwa](https://github.com/lh3/bwa): Necessary for mapping reads onto contigs

3. [bam-readcount](https://github.com/genome/bam-readcount): Used to get per sample base frequencies at each position

4. [samtools] (http://www.htslib.org/download/): Utilities for processing mapped files

5. [CONCOCT](https://github.com/BinPro/CONCOCT): Our own contig binning algorithm

6. [prodigal] (https://github.com/hyattpd/prodigal/releases/): Used for calling genes on contigs

7. [gnu parallel] (http://www.gnu.org/software/parallel/): Used for parallelising rps-blast

8. [standalone blast] (http://www.ncbi.nlm.nih.gov/books/NBK52640/): Need rps-blast

9. COG RPS database: ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/ Cog databases

10. [GFF python parser] (https://github.com/chapmanb/bcbb/tree/master/gff)

and megahit:

package manager:
```
sudo apt-get install bwa samtools

```

<a name="simple_example"/>
##Simple example using Snakemake

###Setup
```
    source activate desman_snakemake
    conda install snakemake
    cd DESMAN_DIR/complete_example_snakemake

    mkdir data; cd data    
    wget https://www.dropbox.com/s/l6g3culvibym8g7/Example.tar.gz
    tar -xvzf Example.tar.gz
    cd ..
```

Make sure you edit the 'config.json' file to suit your system:

```
    less config.json
```

To show what rules will be executed in order to create a certain file, run:

```
    snakemake --dryrun Concoct/Output/clustering_gt1000.csv
```

### Run up to CONCOCT confusion matrix plot

In this step, we'll assemble the reads, run concoct including all the preprocessing of the data and taxonomically assign the contigs.

```
    snakemake --cores 16 confusion_matrix_plot --resources memory=128000000000
```

### Extract clusters associated with E. Coli

In this dataset there is several E. Coli strains present, which CONCOCT is not able to resolve properly.
So we'll use the taxonomic assignment of the contigs to select clusters to use as input for DESMAN. 

Adapt these commands to the clusters that you see are associated with E. Coli in your taxassign_conf.pdf.
 
```
    mkdir Split
    cd Split
    $DESMAN/scripts/SplitClusters.pl ../Contigs/final_contigs_c10K.fa ../Concoct/Output/clustering_gt1000.csv
    cat Cluster14/Cluster14.fa Cluster3/Cluster3.fa Cluster5/Cluster5.fa Cluster19/Cluster19.fa > ../Contigs/ClusterEC.fa 
    cd ..
```

### Run Desman on identified E. Coli contigs

To generate a plot showing mean posterior deviance vs strain number run the command

```
    snakemake --cores 16 RunDesman/ClusterEC/Dev.pdf
```


