# DESMAN complete example using Snakemake

##Managing python environments with Anaconda
DESMAN requires some parts of it to be run with python 3 and some with python 2.7, therefore we strongly recommend the use of Anaconda. Anaconda is a tool to isolate your python installation, which allows you to have multiple parallel installations using different versions of different packages, and gives you a very convenient and fast way to install the most common scientific python packages. Anaconda is free of charge and you can download Anaconda [here](https://store.continuum.io/cshop/anaconda/). Installation instructions can be found [here](http://docs.continuum.io/anaconda/install.html).


##Installation

After installing Anaconda, create an new environment that will contain the python 2.7 installation called 'desman_2_env':

```
    conda create -n desman_2_env python=2.7
```

and one environment that will contain the python 3.5 installation called 'desman_snakemake':


```
    conda create -n desman_snakemake python=3.5
```

In the desman_2_env environment we'll install the concoct dependencies:

```
    source activate desman_2_env
    conda install cython numpy scipy biopython pandas pip scikit-learn
```

and the concoct software:

```
    cd CONCOCT_DIR
    python setup.py install
```

and the desman software:

```
    cd DESMAN_DIR
    python setup.py install
```

To install snakemake from the bioconda channel we'll switch to the desman_snakemake environment:

```
    source activate desman_snakemake
    conda install -c bioconda snakemake
```

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

##Simple example using Snakemake

```
    source activate desman_snakemake
    conda install snakemake
    cd DESMAN_DIR/complete_example_snakemake

    mkdir data; cd data    
    wget https://www.dropbox.com/s/l6g3culvibym8g7/Example.tar.gz
    tar -xvzf Example.tar.gz
    cd ..
    snakemake --list-target-rules
```

This will show a list of possible 'rules' to execute with snakemake.

To show what rules will be executed when running the whole example, run:

```
    snakemake --dryrun classify_contig_nr
```

Make sure you edit the 'config.json' file to suit your system:

```
    less config.json
```

And to execute it run:

```
    snakemake --cores 16 classify_contig_nr --resources memory=128000000000
```

