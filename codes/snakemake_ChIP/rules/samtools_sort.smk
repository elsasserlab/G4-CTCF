rule samtools_sort:
    input:
        "results/{genome}/mapped/{sample}.bam"
    output:
        temp("results/{genome}/mapped/{sample}.sorted.bam")
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:  # Samtools takes additional threads through its option -@
        10     # This value - 1 will be sent to -@.
    conda:
        "../env/conda.yml"
    log:
        "logs/samtools_sort/{genome}/{sample}.log"
    wrapper:
        "0.78.0/bio/samtools/sort"
