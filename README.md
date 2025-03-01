# No evidence that G-quadruplexes promote proximal CTCF binding in vivo

## Figures

Running _plot_figures_ script of the _codes_ directory generates all the figures related to the analysis

## Snakemake workflow for CTCF and G4 Cut&Tag/ChIP-Seq analysis

## Usage

The workflow folder contains the Snakemake environment with the necessary config file and Snakefile for running the Snakemake rules (workflow/rules folder) of the analyses (without normalization). The environment.yaml file listing the necessary packages is recommended for setting up the conda environment (e.g. `conda env create -f environment.yaml` ) before running. The workflow/index/mm10 folder must contain the bowtie2 indexed genome (mm10.*.bt2 formats). Multiple genomes can be provided in the analysis if different fastas are found in the config.yaml file. workflow/raw must contain the fastq pairs. If there are input files, samples and inputs in the workflow/raw folder should be named as "{couple}_chip.fastq.gz" and "{couple}_input.fastq.gz", respectively. Test fastq files are provided in the workflow/raw folder. 

Our ChIP-Seq workflows ran under the versions below: 

    snakemake 7.26
    fastqc 0.12.1
    bowtie2 2.5.1
    samtools 1.17
    multiqc 1.14
    picard 3.0.0
    macs2 2.2.7.1

Python >3.8 and conda (e.g. version 24.1.2) are recommended.

## deepTools figures

The codes folder consists of the deepTools bash scripts used for the deeptools related visualizations. 

# Data

Normalized bigwig files and peaks are stored [here](https://data.mendeley.com/preview/h46bbkndy9?a=cea130ea-7ffc-4484-9841-f9856787ab8b).

# References

[Snakemake documentation](https://snakemake.readthedocs.io/en/stable/)

[deepTools documentation](https://deeptools.readthedocs.io/en/develop/)


