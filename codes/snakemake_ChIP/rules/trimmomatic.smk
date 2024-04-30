rule trimmomatic:
    input:
        "raw/{sample}.fastq.gz"  # input and output can be uncompressed or compressed
    output:
        "results/trimmed/{sample}.trimmed.fastq.gz"
    log:
        "logs/trimmomatic/{sample}.log"
    params:
        # list of trimmers (see manual)
        trimmer=["ILLUMINACLIP:reference_data/NexteraPE-PE.fa.txt:2:30:10"],
        # optional parameters
        extra="",
        # optional compression levels from -0 to -9 and -11
        compression_level="-9"
    threads:
        10
    resources:
        mem_mb=1024
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/trimmomatic/se"
