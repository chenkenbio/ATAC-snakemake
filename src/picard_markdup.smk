
rule picard_markdup_star:
    input:
        outdir / "alignment_star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
    output:
        bam = protected(outdir / "alignment_star/{sample}/{sample}.markdup.bam"),
        bai = protected(outdir / "alignment_star/{sample}/{sample}.markdup.bam.bai"),
        stat = protected(outdir / "alignment_star/{sample}/{sample}.markdup.metrics.txt")
    log:
        outdir / "alignment_star/{sample}/{sample}.markdup.log"
    threads:
        4
    params:
        picard = os.path.expanduser("~/opt/picard.jar")
    shell:
        """
        # add read group
        java -jar {params.picard} AddOrReplaceReadGroups \
            --INPUT {input} \
            --OUTPUT {output.bam}.tmp.bam \
            --RGID {wildcards.sample} \
            --RGLB {wildcards.sample} \
            --RGPL ILLUMINA \
            --RGPU {wildcards.sample} \
            --RGSM {wildcards.sample} &> {log} && \
        java -jar {params.picard} MarkDuplicates \
            --INPUT {output.bam}.tmp.bam \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.stat} &> {log} && \
        samtools index -@ {threads} {output.bam} && \
        rm -f {output.bam}.tmp.bam
        """


rule picard_markdup:
    input:
        outdir / "alignment_hisat2/{sample}/{sample}.sorted.bam"
    output:
        bam = outdir / "alignment_hisat2/{sample}/{sample}.markdup.bam",
        bai = outdir / "alignment_hisat2/{sample}/{sample}.markdup.bam.bai",
        stat = outdir / "alignment_hisat2/{sample}/{sample}.markdup.metrics.txt"
    log:
        outdir / "alignment_hisat2/{sample}/{sample}.markdup.log"
    threads:
        4
    params:
        picard = os.path.expanduser("~/opt/picard.jar")
    shell:
        """
        java -jar {params.picard} MarkDuplicates \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.stat} &> {log} && \
        samtools index -@ {threads} {output.bam} && \
        rm -f {output.bam}.tmp.bam
        """


rule picard_markdup_bowtie2:
    input:
        outdir / "alignment_bowtie2/{sample}/{sample}.sorted.bam"
    output:
        bam = outdir / "alignment_bowtie2/{sample}/{sample}.markdup.bam",
        bai = outdir / "alignment_bowtie2/{sample}/{sample}.markdup.bam.bai",
        stat = outdir / "alignment_bowtie2/{sample}/{sample}.markdup.metrics.txt"
    log:
        outdir / "alignment_bowtie2/{sample}/{sample}.markdup.log"
    threads:
        4
    params:
        picard = os.path.expanduser("~/opt/picard.jar")
    shell:
        """
        java -jar {params.picard} MarkDuplicates \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.stat} &> {log} && \
        samtools index -@ {threads} {output.bam} && \
        rm -f {output.bam}.tmp.bam
        """

rule remove_dup:
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam",
    output:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.rmdup.bam",
        bai = outdir / "alignment_{aligner}/{sample}/{sample}.rmdup.bam.bai"
    threads:
        8
    shell:
        """
        samtools view -F 4 -@ {threads} -F 1024 -b {input.bam} -q 1 > {output.bam} && \
        samtools index -@ {threads} {output.bam}
        """

