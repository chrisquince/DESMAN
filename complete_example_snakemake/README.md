# DESMAN complete example using Snakemake

![alt tag](../desmans.jpg)

## Table of Contents
[Installation](#installation)

[Complete example using Snakemake](#complete_example)

<a name="installation"/>

## Installation

The following was tested on a clean version of Ubuntu 16.04.1 LTS some changes may 
need to be made on other systems.

### Install Dependencies

DESMAN requires some parts of it to be run with python 3 and some with python 2.7,
 therefore we strongly recommend the use of Conda. Conda is a tool to isolate 
 your python installation, which allows you to have multiple parallel installations
  using different versions for different packages, and gives you a very convenient 
  and fast way to install the most common scientific python packages. Additionaly 
  you can use Conda to install all kinds of software into an somewhat isolated 
  environment. Conda is free of charge and you can read more about it [here](https://github.com/conda/conda).

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
    pip install bcbio-gff
```

These items are prerequisities for the installation of DESMAN:

1. *gcc
2. *gsl

The installation procedure varies on different systems,
and described in this README is only how to proceed with a linux (ubuntu) distribution.

A c-compiler, e.g. gcc, is needed to compile the c parts of CONCOCT and DESMAN that use the
GNU Scientific Library GSL. For linux (ubuntu) this is installed through:

```
    sudo apt-get update
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

In the desman_python3_env we will install snakemake, in theory this should be possible within conda 
but we found some conflicts so instead used:

```
source activate desman_python3_env
sudo apt install python3-pip
pip3 install snakemake
```

You will now need to install a number of third party programs. These may be possible through conda 
but mostly we use git or the Ubuntu package managers. We include below our install strategy 
but you may need to adjust these. You will use need to make a local bin directory (mkdir ~/bin) if 
not already present:


1. [megahit](https://github.com/voutcn/megahit): A highly efficient metagenomics assembler currently our default for most studies
    
    ```
    cd ~/repos
    git clone https://github.com/voutcn/megahit
    cd megahit/
    sudo apt-get install zlib1g-dev
    make
    cp megahit* ~/bin
    ```

2. [bwa](https://github.com/lh3/bwa): Necessary for mapping reads onto contigs

    ```
    cd ~/repos
    git clone https://github.com/lh3/bwa.git
    cd bwa; make
    cp bwa ~/bin
    ```

3. [bam-readcount](https://github.com/genome/bam-readcount): Used to get per sample base frequencies at each position

    ```
    cd ~/repos
    sudo apt-get install build-essential git-core cmake zlib1g-dev libncurses-dev patch
    git clone https://github.com/genome/bam-readcount.git
    mkdir bam-readcount-build
    cd bam-readcount-build/
    cmake ../bam-readcount
    make
    cp bin/bam-readcount ~/bin/
    ```

4. [samtools] (http://www.htslib.org/download/): Utilities for processing mapped files. The version 
    available through apt will *NOT* work instead...

    ```
    cd ~/repos
    wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
    tar xvfj samtools-1.3.1.tar.bz2 
    cd samtools-1.3.1/ 
    sudo apt-get install libcurl4-openssl-dev libssl-dev
    ./configure --enable-plugins --enable-libcurl --with-plugin-path=$PWD/htslib-1.3.1
    make all plugins-htslib
    cp samtools ~/bin/  
    ```

5. [bedtools] (http://bedtools.readthedocs.io/en/latest/): Utilities for working with read mappings

    ```
    sudo apt-get install bedtools
    ```

5. [prodigal] (https://github.com/hyattpd/prodigal/releases/): Used for calling genes on contigs

    ```
    wget https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux 
    cp prodigal.linux ~/bin
    chmod +rwx ~/bin/prodigal
    ```

6. [gnu parallel] (http://www.gnu.org/software/parallel/): Used for parallelising rps-blast

    ```
    sudo apt-get install parallel
    ```

7. [standalone blast] (http://www.ncbi.nlm.nih.gov/books/NBK52640/): Need rpsblast from the blast2 package

    ```
    sudo apt-get install blast2
    ```
    
8. [diamond] (https://github.com/bbuchfink/diamond): BLAST compatible accelerated aligner

    ```
    cd ~/repos
    mkdir diamond
    cd diamond
    wget http://github.com/bbuchfink/diamond/releases/download/v0.8.31/diamond-linux64.tar.gz
    tar xzf diamond-linux64.tar.gz
    cp diamond ~/bin/
    ```
    
9. [R] (https://www.r-project.org/) Finally we need R we followed these steps 
[how to install r on linux ubuntu 16-04 xenial xerus](https://www.r-bloggers.com/how-to-install-r-on-linux-ubuntu-16-04-xenial-xerus/):
and installed the additional packages: gplots ggplot2 getopt reshape


We also need some database files versions of which that are compatible with the pipeline 
we have made available through s3. Below we suggest downloading them to a databases directory:

1. COG RPS database: ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/ Cog databases
    
    ```
    mkdir ~/Databases
    cd ~/Databases
    wget https://desmandatabases.s3.climb.ac.uk/rpsblast_cog_db.tar.gz
    tar -xvzf rpsblast_cog_db.tar.gz
    ```
2.  NCBI non-redundant database formatted in old GI format downloaded 02/08/2016 02:07:05. We provide 
this as fasta sequence so that you can diamond format it yourself to avoid any version issues:
    
    ```
    cd ~/Databases
    mkdir NR
    cd NR
    wget https://desmandatabases.s3.climb.ac.uk/nr.faa
    diamond makedb --in nr.faa -d nr
    ```
3. GI to Taxaid and lineage files for the above:
    
    ```
    wget https://desmandatabases.s3.climb.ac.uk/gi_taxid_prot.dmp
    wget https://desmandatabases.s3.climb.ac.uk/all_taxa_lineage_notnone.tsv
    ```




<a name="complete_example"/>
## Complete example using Snakemake

### Setup
```
    source activate desman_python3_env
```

For convenience we will also create environmental variables for the paths to the 
CONCOCT and DESMAN repos:
```
export DESMAN=~/repos/DESMAN/
export CONCOCT=~/repos/CONCOCT/
```

We will work in a new example directory and download the simulated example sequences:

```
mkdir ~/complete_example_snakemake
cd complete_example_snakemake
mkdir data; cd data
wget https://www.dropbox.com/s/l6g3culvibym8g7/Example.tar.gz
tar -xvzf Example.tar.gz
cd ..
```

We now copy in the snakemake config file which will have to be edited to suit your system:
```
cd ~/complete_example_snakemake
cp ~/repos/DESMAN/complete_example_snakemake/config.json .
cp ~/repos/DESMAN/complete_example_snakemake/rules.snake .
cp ~/repos/DESMAN/complete_example_snakemake/Snakefile .
```

Make sure you edit the 'config.json' file to suit your system, it is currently configured based 
on the install above:

```
    less config.json
```

Snakemake works by generating 

To show what rules will be executed in order to create a certain file, run:

```
    snakemake --dryrun Concoct/Output/clustering_gt1000.csv
```

### Run up to CONCOCT confusion matrix plot

In this step, we'll assemble the reads, run concoct including all the preprocessing of the data and taxonomically assign the contigs.

```
    snakemake --cores 8 confusion_matrix_plot --resources memory=6000000
```

### Extract clusters associated with E. Coli

In this dataset there is several E. Coli strains present, which CONCOCT is not able to resolve properly.
So we'll use the taxonomic assignment of the contigs to select clusters to use as input for DESMAN. 

Adapt these commands to the clusters that you see are associated with E. Coli in your taxassign_conf.pdf.
    
    D0   D10   D18    D20     D2
     
```
    mkdir Split
    cd Split
    ~/repos/DESMAN/scripts/SplitClusters.pl ../Contigs/final_contigs_c10K.fa ../Concoct/Output/clustering_gt1000.csv
    cat Cluster0/Cluster0.fa Cluster10/Cluster10.fa Cluster18/Cluster18.fa Cluster20/Cluster20.fa Cluster2/Cluster2.fa > ../Contigs/ClusterEC.fa 
    cd ..
```

### Run Desman on identified E. Coli contigs

To generate a plot showing mean posterior deviance vs strain number run the command

```
    snakemake --cores 8 RunDesman/ClusterEC/Dev.pdf
```


