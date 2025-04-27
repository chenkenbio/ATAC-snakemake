

GENOMES = {
    "hg38": "/home/kenchen/db/gencode/GRCh38/GRCh38.primary_assembly.genome.fa",
    "mm10": "/home/kenchen/db/gencode/GRCm38/GRCm38.primary_assembly.genome.fa",
}
rule blastn_screen:
    wildcard_constraints:
        read="R[12]"
    input:
        outdir / "reads_trimmed/{sample}/{sample}_{read}.trimmed.fastq.gz",
    output:
        fa = temp(outdir / "blastn_screen/{sample}/{sample}_{read}.sample.fa"),
        txt = outdir / "blastn_screen/{sample}/{sample}_{read}.blastn.txt",
    params:
        db = lambda wildcards: GENOMES[SAMPLE_DATA[wildcards.sample]["genome"]],
        num_seqs=10000
    threads:
        8
    log:
        outdir / "blastn_screen/{sample}/{sample}_{read}.blastn.log",
    shell:
        """
        seqtk sample -s 0 {input} {params.num_seqs} | \
            seqtk seq -A > {output.fa} && \
        blastn -query {output.fa} \
            -num_threads {threads} \
            -db {params.db} \
            -outfmt 7 > {output.txt} 2> {log}
        """
        
    

