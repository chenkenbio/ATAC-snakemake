
MACS_GSIZE = {
    "mm10": "mm",
    "hg38": "hs",
}
rule macs_call_narrowpeak_atac:
    input:
        outdir / "alignment_{aligner}/{sample}/{sample}.filtered.bam"
    output:
        npk = outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_peaks.narrowPeak",
        npk_summit = outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_summits.bed",
    params:
        gsize = lambda wildcards: MACS_GSIZE[SAMPLE_DATA[wildcards.sample]["genome"]],
        extra = ' '.join([
            "-f BAMPE",
            "-B",
            "--keep-dup all",
            "--SPMR",
            "-q 0.01",
            # "--nomodel",
            # "--extsize 200",
        ])
    threads:
        1
    shell:
        """
        macs3 callpeak \
            {params.extra} \
            -t {input} \
            -g {params.gsize} \
            -n {wildcards.sample} \
            --outdir {outdir}/macs3_narrowpeak_{wildcards.aligner}/{wildcards.sample} \
            &> {output.npk}.log
        """

rule count_reads_in_peaks:
    wildcard_constraints:
        dupmode = r"rmdup|markdup|filtered"
    input:
        npk = outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_peaks.narrowPeak",
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.{dupmode}.bam"
    output:
        counts = outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_{dupmode}_peaks.counts",
        log = outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_{dupmode}_peaks.counts.log"
    threads:
        16
    shell:
        """
        src/biock_featureCounts.py \
            -a {input.npk} \
            -F BED \
            -s 0 \
            -i {input.bam} \
            -p --countReadPairs \
            -T {threads} \
            -o {output.counts} \
            2> {output.log}
        """

rule get_FRiP:
    input:
        outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_filtered_peaks.counts.log"
    output:
        outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_filtered.FRiP.txt"
    params:
        sample = lambda wildcards: wildcards.sample
    shell:
        """
        res=$(grep 'Successfully assigned alignments :' {input} | awk '{{print $6,$7}}' | sed 's/)//;s/(//;s/ /\t/')
        echo -e "#sample\tFRiP" > {output}
        echo -e "{params.sample}\t$res" >> {output}
        """
