# VirusMap

VirusMap is a RNA-seq analysis tool built using snakemake. I tried using Magic-Blast to look for reads in the SRA, but the connection is not stable for searching large datasets. VirusMap will download reads and map locally. Aligment is performed by SRA with a specifc target genome. I'll probally go back to Magic-Blast for more general local searches.

## Directed acyclic graph aka the pipeline

![virusMAP DAG](./rulegraph.png "virusMAP DAG")

## Setup

First clone this repo...

$ `bash setup.sh`

## Run

Set `SRA_IDS` and `VIRAL_GENBANK_IDS` in Snakefile

$ `snakemake`

## Pipleline

### Downloading data

- Curl is backup
- Fasta-dump freaked out and downloaded 32Gb of data instead of the 4Gb reads...
- Biopython - it's python

### Pre-align-QC

- FASTQC
- Trim Galore!

### Mapping

- STAR
- Maybe Magic-Blast

### Post-Align-QC

- Qualimap

### Reporting

- MultiQC

## More info

The ncbi has a good SnakeMake tutorial <https://nbis-reproducible-research.readthedocs.io/en/latest/snakemake/>

Create conda enviroment without prefix:  
`conda env export | grep -v "^prefix: " > environment.yml`