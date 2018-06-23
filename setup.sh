#!/usr/bin/bash

# add channels
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels r
conda config --add channels bioconda

# create conda env
conda create --yes --name virusMAP 

#install dependancies
echo('Installing dependancies, this will take a while...')
echo('I turned on debug so you have something to watch ;)')

sleep 3

conda install --yes --debug -n virusMAP biopython=1.70 
conda install --yes --debug -n virusMAP fastqc=0.11.7
conda install --yes --debug -n virusMAP snakemake=5.1.4
conda install --yes --debug -n virusMAP trim-galore=0.4.5
conda install --yes --debug -n virusMAP star=2.6.0c 
conda install --yes --debug -n virusMAP qualimap=2.2.2a
conda install --yes --debug -n virusMAP multiqc=1.5

echo('To activate the virusMAP environment, use')
echo(' $ conda activate virusMAP')