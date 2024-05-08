rule filtering_unmapped:
    input:
        "results/{genome}/mapped/{sample}.sorted.bam"
    output:
        temp("results/{genome}/mapped/{sample}.sorted.filtered.bam")
    log:
        "logs/filtering_unmapped/{genome}/{sample}.log"
    params:
        extra="-F 0x4" # optional params string
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/samtools/view"