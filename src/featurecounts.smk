
for name, settings in config["featureCounts"].items():
    if "saf" in settings:
        anno_fn = settings["saf"]
        feature = "-F SAF"
    else:
        anno_fn = settings["gtf"]
        feature = ' '.join([
            "-F GTF",
            '-t', settings["feature"],
            '-g', settings["metafeature"]
        ])
    rule:
        wildcard_constraints:
            dupmode = "markdup|rmdup|filtered"
        priority: 1000
        name: f"featureCounts_{name}",
        input:
            bam = outdir / "alignment_{aligner}/{sample}/{sample}.{dupmode}.bam",
        output:
            counts = outdir / "featureCounts_{aligner}" / f"{name}.{{dupmode}}" / "{sample}/{sample}.counts",
        params:
            anno = anno_fn,
            pe = "-p --countReadPairs" if settings["paired_end"] else "",
            feature = feature,
            strand = settings["strand"], # 0: unstranded, 1: stranded, 2: reverse stranded
            extra = ' '.join(settings.get("extra", [])),
            name = lambda wildcards: wildcards.sample,
        threads:
            8
        shell:
            """
            src/biock_featureCounts.py \
                -i {input.bam} \
                -T {threads} \
                -a {params.anno} \
                -o {output.counts} \
                -s {params.strand} \
                --name {params.name} \
                {params.feature} \
                {params.pe} \
                {params.extra} &> {output.counts}.log
            """

rule merge_featureCounts:
    wildcard_constraints:
        dupmode = "markdup|rmdup|filtered"
    input:
        expand(outdir / "featureCounts_{{aligner}}/{{name}}.{{dupmode}}/{sample}/{sample}.counts", sample=SAMPLE_DATA.keys())
    output:
        outdir / "featureCounts_{aligner}" / "{name}.{dupmode}" / "{name}_merged.counts"
    params:
        labels = list(SAMPLE_DATA.keys())
    shell:
        """
        src/biock_merge-featureCounts.py \
            -i {input} \
            -l {params.labels} \
            -p {outdir}/featureCounts_{wildcards.aligner}/{wildcards.name}.{wildcards.dupmode}/{wildcards.name}_merged
        """
 
    
# usage: biock_merge-featureCounts.py [-h] -i FEATURECOUNTS [FEATURECOUNTS ...]
#                                     [-l LABEL [LABEL ...]] -p PREFIX
#                                     [--ignore-summary] [--no-header]
#                                     [--verbose]
# 
# options:
#   -h, --help            show this help message and exit
#   -i FEATURECOUNTS [FEATURECOUNTS ...], --featureCounts FEATURECOUNTS [FEATURECOUNTS ...]
#                         featureCounts count file
#   -l LABEL [LABEL ...], --label LABEL [LABEL ...]
#                         labels for each featureCounts file
#   -p PREFIX, --prefix PREFIX
#                         prefix for output file
#   --ignore-summary      Ignore summary file
#   --no-header           ignore header in featureCounts file
#   --verbose             Print verbose output
    
       

