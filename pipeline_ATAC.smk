#!/usr/bin/env snakemake

## Author: Ken Chen <chenkenbio>

import os
import sys
import yaml
import json
from pathlib import Path

include: "src/helper.smk"

configfile: "demo/config.yaml"

outdir = Path(config["outdir"])
SAMPLE_DATA = load_sample_info(config["samples"])

assert len(config) > 0

config["featureCounts"] = yaml.load(open(config["featureCounts_config"]), Loader=yaml.FullLoader)

rule all:
    input:
        expand(outdir / "reads_raw/{sample}/{sample}_{read}_fastqc.html", sample=SAMPLE_DATA.keys(), read=["R1", "R2"]),
        expand(outdir / "reads_trimmed/{sample}/{sample}_{read}.trimmed.fastq.gz", sample=SAMPLE_DATA.keys(), read=["R1", "R2"]),
        expand(outdir / "reads_trimmed/{sample}/{sample}_{read}.trimmed_fastqc.html", sample=SAMPLE_DATA.keys(), read=["R1", "R2"]),
        expand(outdir / "fastq_screen/{sample}/{sample}_{read}.trimmed_screen.txt", sample=SAMPLE_DATA.keys(), read=["R1", "R2"]),
        expand(outdir / "blastn_screen/{sample}/{sample}_{read}.blastn.txt", sample=SAMPLE_DATA.keys(), read=["R1", "R2"]),
        expand(outdir / "alignment_bowtie2/{sample}/{sample}.filtered.bam.bai", sample=SAMPLE_DATA.keys()),
        expand(
            outdir / "alignment_bowtie2/{sample}/{sample}.{metric}.txt",
            sample=SAMPLE_DATA.keys(),
            metric=["samtools_stats", "InsertSizeMetrics", "LibraryComplexity"],
        ),
        expand(
            outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_filtered_peaks.counts", 
            aligner=["bowtie2"], 
            sample=SAMPLE_DATA.keys(),
        ),
        expand(outdir / "macs3_narrowpeak_{aligner}/{sample}/{sample}_filtered.FRiP.txt", aligner=["bowtie2"], sample=SAMPLE_DATA.keys()),
        expand(
            outdir / "featureCounts_bowtie2/{name}.filtered/{sample}/{sample}.counts",
            sample=SAMPLE_DATA.keys(),
            name=config["featureCounts"].keys(),
        ),
        expand(
            outdir / "featureCounts_bowtie2/{name}.filtered/{name}_merged.counts",
            name=config["featureCounts"].keys(),
        ),
        expand(
            outdir / "coverage_{aligner}_{norm}_mapq{mapq}" / "{sample}/{sample}.filtered.mapq{mapq}.{norm}.bw",
            aligner=["bowtie2"],
            norm=["None", "RPKM"],
            mapq=[1],
            sample=SAMPLE_DATA.keys(),
        ),
        expand(
            outdir / "coverage_{aligner}_{norm}_mapq{mapq}" / "{sample}/{sample}.filtered.mapq{mapq}.{norm}.gene-matrix.gz",
            aligner=["bowtie2"],
            norm=["None", "RPKM"],
            mapq=[1],
            sample=SAMPLE_DATA.keys(),
        )
        
include: "src/link_fastq.smk"
include: "src/fastp_v0.smk"
include: "src/fastq_screen.smk"
include: "src/blastn_screen.smk"
include: "src/fastqc.smk"
include: "src/atac_bowtie2.smk"
include: "src/filter_bam.smk"
include: "src/picard_markdup.smk"
include: "src/macs3_atac.smk"
include: "src/alignment_qc_v0.smk"
include: "src/bamcoverage_atac.smk"
include: "src/featurecounts.smk"
