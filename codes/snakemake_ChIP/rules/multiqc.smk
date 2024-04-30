rule multiqc:
    input:
        expand("results/{genome}/samtools_flagstat/{sample}.sorted.filtered.marked.flagstat", sample=samples, genome=genomes)
    output:
        "results/qc/multiqc.html"
    params:
        ""  # Optional: extra parameters for multiqc.
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/multiqc"