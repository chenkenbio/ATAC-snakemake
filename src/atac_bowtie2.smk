rule bowtie2_atac_pe:
    input:
        R1 = outdir / "reads_trimmed/{sample}/{sample}_R1.trimmed.fastq.gz",
        R2 = outdir / "reads_trimmed/{sample}/{sample}_R2.trimmed.fastq.gz"
    output:
        unsorted = temp(outdir / "alignment_bowtie2/{sample}/{sample}.bam"),
        bam = temp(outdir / "alignment_bowtie2/{sample}/{sample}.sorted.bam"),
        bai = temp(outdir / "alignment_bowtie2/{sample}/{sample}.sorted.bam.bai"),
    params:
        index = lambda wildcards: config["bowtie2_index"][SAMPLE_DATA[wildcards.sample]["genome"]],
        read_group = lambda wildcards: f"--rg-id {wildcards.sample} --rg SM:{wildcards.sample}",
        extra = ' '.join([
            "--very-sensitive",
        ])
    threads:
        24
    log:
        out = outdir / "alignment_bowtie2/{sample}/{sample}.log"
    shell:
        """
        bowtie2 \
            -p {threads} \
            -x {params.index} \
            {params.read_group} \
            {params.extra} \
            -X 2000 \
            --un-conc-gz {output.unsorted}.unmapped.gz \
            -1 {input.R1} -2 {input.R2} 2> {log.out} | samtools view -b -@ {threads} > {output.unsorted} && \
        samtools sort -@ {threads} -o {output.bam} {output.unsorted} && \
        samtools index -@ {threads} {output.bam}
        """
        

