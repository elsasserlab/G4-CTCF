rule fastqc_trim:
    input:
        "results/trimmed/{sample}.trimmed.fastq.gz"
    output:
        html="results/qc/fastqc/trimmed/{sample}.html",
        zip="results/qc/fastqc/trimmed/{sample}_fastqc.zip" # the suffix _fastqc.zip is necessary for multiqc to find the file. If not using multiqc, you are free to choose an arbitrary filename
    params: "--quiet"
    log:
        "logs/fastqc/trimmed/{sample}.log"
    threads: 10
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/fastqc"
