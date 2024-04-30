rule fastqc:
    input:
        "raw/{sample}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params:
    log:
        "logs/fastqc/{sample}.log"
    threads: 5
    conda:
        "../env/conda.yml"
    wrapper:
	    "v1.26.0/bio/fastqc"
