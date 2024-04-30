rule mark_duplicates:
    input:
        sample="results/{genome}/mapped/{sample}.sorted.filtered.bam"
    # optional to specify a list of BAMs; this has the same effect
    # of marking duplicates on separate read groups for a sample
    # and then merging
    output:
        bam="results/{genome}/dedup/{sample}.sorted.filtered.marked.bam",
        metrics="results/{genome}/dedup/{sample}.metrics.txt"
    log:
        "logs/dedup/{genome}/{sample}.log"
    params:
        extra="REMOVE_DUPLICATES=true"
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/picard/markduplicates"