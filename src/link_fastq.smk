ruleorder: link_fastq_pe > link_fastq_se

rule link_fastq_pe:
    wildcard_constraints:
        sample = '|'.join([sn for sn in SAMPLE_DATA.keys() if 'R2' in SAMPLE_DATA[sn]])
    input:
        r1 = lambda wildcards: SAMPLE_DATA[wildcards.sample]['R1'],
        r2 = lambda wildcards: SAMPLE_DATA[wildcards.sample].get('R2', '/')
    output:
        r1 = outdir / "reads_raw/{sample}/{sample}_R1.fastq.gz",
        r2 = outdir / "reads_raw/{sample}/{sample}_R2.fastq.gz"
    threads:
        1
    shell:
        "ln -sf {input.r1} {output.r1} && ln -sf {input.r2} {output.r2}"

rule link_fastq_se:
    wildcard_constraints:
        sample = '|'.join([sn for sn in SAMPLE_DATA.keys() if ('R2' not in SAMPLE_DATA[sn] or SAMPLE_DATA[sn]['R2'] is None)])
    input:
        r1 = lambda wildcards: SAMPLE_DATA[wildcards.sample]['R1'],
    output:
        r1 = outdir / "reads_raw/{sample}/{sample}_R1.fastq.gz",
    threads:
        1
    shell:
        "ln -sf {input.r1} {output.r1}"
