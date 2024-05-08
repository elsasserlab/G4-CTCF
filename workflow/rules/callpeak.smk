rule callpeak:
    input:
        bai="results/{genome}/dedup/{couple}_input.sorted.filtered.marked.bam.bai",
        treatment="results/{genome}/dedup/{couple}_chip.sorted.filtered.marked.bam",   # required: treatment sample(s)
        control="results/{genome}/dedup/{couple}_input.sorted.filtered.marked.bam"      # optional: control sample(s)
    output:
        # all output-files must share the same basename and only differ by it's extension
        # Usable extensions (and which tools they implicitly call) are listed here:
        #         https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/macs2/callpeak.html.
        multiext("results/{genome}/callpeak/",
                 "{couple}_peaks.xls",   ### required
                 ### optional output files
                 "{couple}_peaks.narrowPeak",
                 "{couple}_summits.bed"
                 )
    log:
        "logs/macs2/{genome}/callpeak_{couple}.log"
    params:
        "-f BAM -g hs --nomodel -q 0.01"
    conda:
        "../env/conda.yml"
    wrapper:
        "0.78.0/bio/macs2/callpeak"
