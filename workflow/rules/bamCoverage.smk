rule bamCoverage:
    input:
        bam="results/{genome}/dedup/{sample}.sorted.filtered.marked.bam",
        bai="results/{genome}/dedup/{sample}.sorted.filtered.marked.bam.bai"
    output:
        "results/bamCoverage/{genome}/{sample}.bw"
    params:
        lambda wildcards:"--effectiveGenomeSize " + config["effectiveGenomeSize"]+ " --blackListFileName " + config["blacklist_file"]
    log:
        "logs/bamCoverage/{genome}/{sample}.log"
    threads: 10
    conda:
        "../env/conda.yml"
    shell:
        "bamCoverage -b {input.bam} -o {output} -of bigwig -p {threads} {params} 2> {log}"

