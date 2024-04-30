rule samtools_stats:
    input:
        "results/{genome}/dedup/{sample}.sorted.filtered.marked.bam"
    output:
        "results/{genome}/samtools_stats/{sample}.txt"
    params:
        extra=""                       # Optional: extra arguments.
    log:
        "logs/samtools_stats/{genome}/{sample}.log"
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/samtools/stats"