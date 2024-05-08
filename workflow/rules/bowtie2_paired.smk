rule bowtie2:
    input:
        sample=["raw/{sample}_1.fastq.gz", "raw/{sample}_2.fastq.gz"],
        index=lambda wildcards:multiext(
        "index/"+ wildcards.genome+"/"+wildcards.genome,
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    output:
        temp("results/{genome}/mapped/{sample}.bam")
    log:
        "logs/bowtie2/{genome}/{sample}.log"
    params:
        index=lambda wildcards:"index/"+ wildcards.genome+"/"+wildcards.genome,
        extra=""  # optional parameters
    threads: 20  # Use at least two threads
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/bowtie2/align"
