rule callpeak_no_couples:
    input:
        bai="results/{genome}/dedup/{sample}.sorted.filtered.marked.bam.bai",
        treatment="results/{genome}/dedup/{sample}.sorted.filtered.marked.bam"   # required: treatment sample(s)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/{genome}/callpeak/",
                 "{sample}_peaks.xls",   ### required
                 ### optional output files
                 "{sample}_peaks.narrowPeak",
                 "{sample}_summits.bed"
                 )
    log:
        "logs/macs2/{genome}/callpeak_{sample}.log"
    params:
        "-f BAM -g hs --nomodel -q 0.01"
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/macs2/callpeak"
