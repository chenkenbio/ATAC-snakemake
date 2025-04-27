#!/bin/bash

# for min_mapq in config["min_mapq"]:
#     for dupmode in ["markdup", "rmdup"]:
#         rule: 
#             name: f"filter_{dupmode}_q{min_mapq}"
#             input:
#                 outdir / "filtered_alignment_{aligner}/{sample}"/ ("{sample}" + f"{dupmode}.bam")
#             output:
#                 bam = outdir / "filtered_alignment_{aligner}_q{minq}/{sample}" / ("{sample}" + f".{dupmode}.q{min_mapq}.bam"),
#                 bai = outdir / "filtered_alignment_{aligner}_q{minq}/{sample}" / ("{sample}" + f".{dupmode}.q{min_mapq}.bam.bai")
#             params:
#                 min_mapq = min_mapq
#             threads:
#                 4
#             shell:
#                 """
#                 samtools view -@ {threads} -q {params.min_mapq} {input} -b > {output.bam} && \
#                 samtools index -@ {threads} {output.bam}
#                 """
BLACKLIST_DICT = {
    "hg38": "/home/kenchen/db/blacklist/hg38-blacklist.v2.bed",
    "mm10": "/home/kenchen/db/blacklist/mm10-blacklist.v2.bed"
}

rule filter_bam:
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam",
    output:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.filtered.bam",
        bai = outdir / "alignment_{aligner}/{sample}/{sample}.filtered.bam.bai"
    params:
        bl = lambda wildcards: BLACKLIST_DICT[SAMPLE_DATA[wildcards.sample]["genome"]],
    threads:
        8
    shell:
        """
        src/biock_filter_bam.py \
            -i {input.bam} \
            -o {output.bam} \
            -t {threads} \
            --rmdup \
            --blacklist {params.bl} \
            --exclude-chromosomes chrM chrMt Mt && \
        samtools index -@ {threads} {output.bam}
        """
#usage: biock_filter_bam.py [-h] [--rmdup] [--min-mapq MIN_MAPQ] [-t THREADS]
#                           [--blacklist BLACKLIST]
#                           [--exclude-chromosomes EXCLUDE_CHROMOSOMES [EXCLUDE_CHROMOSOMES ...]]
#                           [--chr-only]
#                           bam
#biock_filter_bam.py: error: the following arguments are required: bam
        



# rule filter_markdup:
#     input:
#         outdir / "alignment_{aligner}/{sample}/{sample}.markdup.bam"
#     output:
#         bam = outdir / "alignment_{aligner}_q{minq}/{sample}" / ("{sample}" + f".markdup.q{min_mapq}.bam"),
#         bai = outdir / "alignment_{aligner}_q{minq}/{sample}" / ("{sample}" + f".markdup.q{min_mapq}.bam.bai")
#     params:
#         min_mapq = min_mapq
#     threads:
#         4
#     shell:
#         """
#         samtools view -@ {threads} -q {params.min_mapq} {input} -b > {output} && \
#         samtools index -@ {threads} {output}
#         """
# 
# 
# 
