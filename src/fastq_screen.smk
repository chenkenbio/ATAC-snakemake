import os

rule fastq_screen:
    wildcard_constraints:
        read = "R1|R2",
    input:
        outdir / "reads_trimmed/{sample}/{sample}_{read}.trimmed.fastq.gz",
    output:
        outdir / "fastq_screen/{sample}/{sample}_{read}.trimmed_screen.txt",
    params:
        bin = config.get("fastq_screen", {}).get("bin", os.path.expanduser("~/opt/FastQ-Screen-0.16.0/fastq_screen")),
        config = config.get("fastq_screen", {}).get("config", os.path.expanduser("~/opt/FastQ-Screen-0.16.0/fastq_screen.conf")),
    log:
        outdir / "fastq_screen/{sample}_{read}.log",
    threads: 
        4,
    shell:
        """
        rm -f {outdir}/fastq_screen/{wildcards.sample}/{wildcards.sample}_{wildcards.read}* 
        {params.bin} \
            --aligner bowtie2 \
            --subset 200000 \
            --threads {threads} \
            --conf {params.config} \
            --outdir {outdir}/fastq_screen/{wildcards.sample} {input} &> {log}
        """
