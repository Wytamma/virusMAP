#!/usr/bin/bash

# create conda env
conda create --yes --name virusMAP 

# activate conda env
source activate virusMAP

# add channels
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

#install dependancies
conda install biopython=1.70 trim-galore=0.4.5 snakemake=5.1.4 star=2.6.0c fastqc=0.11.7 qualimap=2.2.2a multiqc=1.5