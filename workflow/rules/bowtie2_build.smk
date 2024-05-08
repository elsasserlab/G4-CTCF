rule bowtie2_build:
    input:
        reference = lambda wildcards:config["genome"][wildcards.genome]
    output:
        multiext(
        "index/{genome}/{genome}",
        ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
    log:
        "logs/bowtie2_build/{genome}.log"
    params:
        extra=""  # optional parameters
    threads: 40
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/bowtie2/build"
