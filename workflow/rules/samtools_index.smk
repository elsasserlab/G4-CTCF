rule samtools_index:
    input:
        "results/{genome}/dedup/{sample}.sorted.filtered.marked.bam"
    output:
        "results/{genome}/dedup/{sample}.sorted.filtered.marked.bam.bai"
    log:
        "logs/samtools_index/{genome}/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        10     # This value - 1 will be sent to -@
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/samtools/index"
