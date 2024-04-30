rule samtools_flagstat:
    input:
        "results/{genome}/dedup/{sample}.sorted.filtered.marked.bam"
    output:
        "results/{genome}/samtools_flagstat/{sample}.sorted.filtered.marked.flagstat"
    conda:
        "../env/conda.yml"
    log:
        "logs/samtools_flagstat/{genome}/{sample}.log"
    wrapper:
        "0.78.0/bio/samtools/flagstat"