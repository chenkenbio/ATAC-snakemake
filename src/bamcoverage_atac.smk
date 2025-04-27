GENOME_SIZES = {
    "hg38": 2913022398,
    "mm10": 2652783500
}
BLACKLIST = {
    "hg38": "/home/kenchen/db/blacklist/hg38-blacklist.v2.bed",
    "mm10": "/home/kenchen/db/blacklist/mm10-blacklist.v2.bed"
}

rule bamcoverage_dnaseq:
    wildcard_constraints:
        dupmode = "markdup|rmdup|filtered",
        norm = "BPM|None|CPM|RPGC|RPKM"
    input:
        bam = outdir / "alignment_{aligner}/{sample}/{sample}.{dupmode}.bam",
    output:
        outdir / "coverage_{aligner}_{norm}_mapq{mapq}" / "{sample}/{sample}.{dupmode}.mapq{mapq}.{norm}.bw"
    log:
        temp(outdir / "coverage_{aligner}_{norm}_mapq{mapq}" / "{sample}/{sample}.{dupmode}.mapq{mapq}.{norm}.log")
    params:
        bs = config.get("bamCoverage", {}).get("binSize", 50),
        effsize = lambda wildcards: GENOME_SIZES[SAMPLE_DATA[wildcards.sample]["genome"]],
        extra = lambda wildcards: ' '.join([
            "--ignoreForNormalization chrX chrY chrM",
            f"--normalizeUsing {wildcards.norm}",
            f"--minMappingQuality {wildcards.mapq}",
        ]),
        blacklist = lambda wildcards: BLACKLIST[SAMPLE_DATA[wildcards.sample]["genome"]],
    threads:
        8
    shell:
        """
        bamCoverage \
            -b {input.bam} \
            -o {output} \
            -p {threads} \
            --binSize {params.bs} \
            --blackListFileName {params.blacklist} \
            --effectiveGenomeSize {params.effsize} \
            {params.extra} &> {log}
        """

GTF = {
    "hg38": "/home/kenchen/db/gencode/GRCh38/v46/gencode.v46.annotation.canonical.gtf",
    "mm10": "/home/kenchen/db/gencode/GRCm38/vM25/gencode.vM25.annotation.gtf"
}

rule computematrix_per_sample:
    priority: 100
    input:
        bw = outdir / "coverage_{aligner}_{norm}_mapq{mapq}/{sample}/{sample}.{dupmode}.mapq{mapq}.{norm}.bw"
    output:
        outdir / "coverage_{aligner}_{norm}_mapq{mapq}/{sample}/{sample}.{dupmode}.mapq{mapq}.{norm}.gene-matrix.gz"
    params:
        gtf = lambda wildcards: GTF[SAMPLE_DATA[wildcards.sample]["genome"]],
        labels = lambda wildcards: f"{wildcards.sample}",
        bl = lambda wildcards: BLACKLIST[SAMPLE_DATA[wildcards.sample]["genome"]]
    log:
        temp(outdir / "coverage_{aligner}_{norm}_mapq{mapq}/{sample}/{sample}.{dupmode}.mapq{mapq}.{norm}.gene-matrix.log")
    threads:
        16
    shell:
        """
        computeMatrix scale-regions \
            -S {input.bw} \
            -R {params.gtf} \
            --samplesLabel {params.labels} \
            -a 3000 \
            -b 3000 \
            --regionBodyLength 6000 \
            --skipZeros \
            -o {output} \
            --blackListFileName {params.bl} \
            --missingDataAsZero \
            --binSize 100 \
            --numberOfProcessors {threads} &> {log}
        """



def collate_bams(wildcards):
    bams = list()
    for sn, data in SAMPLE_DATA.items():
        if data["cell"] == wildcards.cell:
            bams.append(outdir / f"coverage_{wildcards.aligner}_{wildcards.norm}" / f"{sn}/{sn}.bw")
    return bams


rule computematrix:
    priority: 100
    input:
        bam = collate_bams
    output:
        outdir / "coverage_{aligner}_{norm}/{cell}.matrix.mat.gz"
    params:
        gtf = "/home/kenchen/db/gencode/GRCh38/v46/gencode.v46.annotation.gtf",
        labels = lambda wildcards: ' '.join([x for x, data in SAMPLE_DATA.items() if data["cell"] == wildcards.cell]),
        bl = "/home/kenchen/db/blacklist/hg38-blacklist.v2.bed"
    log:
        temp(outdir / "coverage_{aligner}_{norm}/{cell}.matrix.log")
    threads:
        32
    shell:
        """
        computeMatrix scale-regions \
            -S {input.bam} \
            -R {params.gtf} \
            --samplesLabel {params.labels} \
            -a 3000 \
            -b 3000 \
            --regionBodyLength 6000 \
            --skipZeros \
            -o {output} \
            --blackListFileName {params.bl} \
            --missingDataAsZero \
            --binSize 100 \
            --numberOfProcessors {threads} &> {log}
        """

rule multibigwigsummary:
    input:
        bam = collate_bams
    output:
        outdir / "coverage_{aligner}_{norm}/{cell}.summary.npz"
    params:
        labels = lambda wildcards: ' '.join([x for x, data in SAMPLE_DATA.items() if data["cell"] == wildcards.cell])
    log:
        temp(outdir / "coverage_{aligner}_{norm}/{cell}.summary.log")
    threads:
        32
    shell:
        """
        multiBigwigSummary bins \
            -b {input.bam} \
            --labels {params.labels} \
            -o {output} \
            --binSize 5000 \
            --numberOfProcessors {threads} &> {log}
        """

