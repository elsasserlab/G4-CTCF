bw files for now: https://filesender.renater.fr/?s=download&token=719c78b5-f66f-4a14-a136-05924840cd71

parameters used: 

computeMatrix
binSize: 10
bef: 2000
blacklist: mm10-blacklist.v2.bed

plotHeatmap
colorMap: YlOrRd


For the ChIP analyses, use snakemake_ChIP

Generic snakemake workflow for ChIP-seq analysis (no normalization)

Dependencies:

    snakemake (e.g. 7.26)
    fastqc (e.g. 0.12.1)
    bowtie2 (e.g. 2.5.1)
    samtools (e.g. 1.17)
    multiqc (e.g. 1.14)
    picard (e.g. 3.0.0)
    macs2 (e.g. 2.2.7.1)

There's a config file for the genome name and fasta if bowtie2 index are not provided, as well as for deeptools bamCoverage parameters.

Naming rules: Sample FASTQs should be put in the raw folder, gzipped. They should be named "raw/{sample].fastq.gz"

If there are inputs, ChIP samples should be named "raw/{couple}_chip.fastq.gz". Input samples should be named "raw/{couple}_input.fastq.gz".

Multiple genomes can be provided in the analysis if different fastas are found in the config.yaml file.
