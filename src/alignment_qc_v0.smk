import os



RSEQC_DB_PATH = os.path.expanduser("~/db/rseqc")
RSEQC_DB = {
    "danRer11": f"{RSEQC_DB_PATH}/danRer11_refseq.bed",
    "hg38": f"{RSEQC_DB_PATH}/hg38_refSeq.bed",
    "mm10": f"{RSEQC_DB_PATH}/mm10_RefSeq.bed",
    "mm39": f"{RSEQC_DB_PATH}/mm39_refseq.bed",
    "TAIR10": f"{RSEQC_DB_PATH}/TAIR10.refseq.bed",
}
rule readdistribution:
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
    output:
        outdir / "alignment_{aligner}/{sample}/{sample}.read_distribution.txt"
    log:
        outdir / "alignment_{aligner}/{sample}/{sample}.read_distribution.log"
    params:
        read_distribution = os.path.expanduser("~/local/conda/rseqc/bin/read_distribution.py"),
        ref = lambda wildcards: RSEQC_DB[SAMPLE_DATA[wildcards.sample]["genome"]]
    threads:
        1
    shell:
        """
        {params.read_distribution} -i {input.bam} -r {params.ref} > {output} 2> {log}
        """

rule infer_experiment:
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
    output:
        outdir / "alignment_{aligner}/{sample}/{sample}.infer_experiment.txt"
    log:
        outdir / "alignment_{aligner}/{sample}/{sample}.infer_experiment.log"
    params:
        infer_experiment = os.path.expanduser("~/local/conda/rseqc/bin/infer_experiment.py"),
        ref = lambda wildcards: RSEQC_DB[SAMPLE_DATA[wildcards.sample]["genome"]],
    threads:
        1
    shell:
        """
        {params.infer_experiment} -i {input.bam} -r {params.ref} > {output} 2> {log}
        """

rule samtools_stats:
    priority: 10
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
    output:
        outdir / "alignment_{aligner}/{sample}/{sample}.samtools_stats.txt"
    threads:
        2
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output}
        """

REF_FLAT = {
    "mm10": os.path.expanduser("~/db/picard/mm10_refFlat.txt"),
    "hg38": os.path.expanduser("~/db/picard/hg38_refFlat.txt"),
    "TAIR10": "",
}
#PICARD_RRNA_LIST = {
#    "mm10": os.path.expanduser("~/db/picard/GRCm38-vM25_rRNA.list"),
#    "hg38": os.path.expanduser("~/db/picard/GRCh38-v46_rRNA.list"),
#    "TAIR10": "",
#}
#rule picard_rnaseqmetrics:
#    priority: 50
#    input:
#        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
#    output:
#        outdir / "alignment_{aligner}/{sample}/{sample}.RnaSeqMetrics.txt"
#    log:
#        temp(outdir / "alignment_{aligner}/{sample}/{sample}.RnaSeqMetrics.log")
#    params:
#        picard = os.path.expanduser("~/opt/picard.jar"),
#        ref_flat = lambda wildcards: REF_FLAT[SAMPLE_DATA[wildcards.sample]["genome"]],
#        ribosomal_intervals = lambda wildcards: PICARD_RRNA_LIST[SAMPLE_DATA[wildcards.sample]["genome"]]
#    threads:
#        4
#    shell:
#        """
#        biock_picard-collectrnaseqmetrics.py \
#            -i {input.bam} \
#            -o {output} \
#            -r {params.ref_flat} \
#            -s SECOND_READ_TRANSCRIPTION_STRAND \
#            -rRNA {params.ribosomal_intervals} > {log} 2>&1
#        """

rule picard_collectinsertsize:
    input:
        outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
    output:
        outdir / "alignment_{aligner}/{sample}/{sample}.InsertSizeMetrics.txt"
    log:
        temp(outdir / "alignment_{aligner}/{sample}/{sample}.InsertSizeMetrics.log")
    params:
        picard = os.path.expanduser("~/opt/picard.jar")
    threads:
        4
    shell:
        """
        java -jar {params.picard} CollectInsertSizeMetrics \
            I={input} \
            O={output} \
            H={outdir}/alignment_{wildcards.aligner}/{wildcards.sample}/{wildcards.sample}.InsertSizeMetrics.pdf \
            &> {log}
        """
    

rule picard_estimatelibrarycomplexity:
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
    output:
        outdir / "alignment_{aligner}/{sample}/{sample}.LibraryComplexity.txt"
    log:
        temp(outdir / "alignment_{aligner}/{sample}/{sample}.LibraryComplexity.log")
    params:
        picard = os.path.expanduser("~/opt/picard.jar")
    threads:
        4
    shell:
        """
        java -jar {params.picard} EstimateLibraryComplexity \
            I={input.bam} \
            O={output} \
            2> {log}
        """
