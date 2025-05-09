ruleorder: fastqc_posttrim > fastqc_posttrim_se
ruleorder: fastqc_pretrim > fastqc_pretrim_se

rule fastqc_pretrim:
    input:
        r1 = outdir / "reads_raw/{sample}/{sample}_R1.fastq.gz",
        r2 = outdir / "reads_raw/{sample}/{sample}_R2.fastq.gz"
    output:
        r1 = outdir / "reads_raw/{sample}/{sample}_R1_fastqc.html",
        r2 = outdir / "reads_raw/{sample}/{sample}_R2_fastqc.html"
    log:
        temp(outdir / "reads_raw/{sample}.fastqc.log")
    threads:
        2
    shell:
        """
            fastqc \
                -t {threads} \
                -o "{outdir}/reads_raw/{wildcards.sample}" \
                {input.r1} {input.r2} &> {log}
        """


rule fastqc_posttrim:
    input:
        r1 = outdir / "reads_trimmed/{sample}/{sample}_R1.trimmed.fastq.gz",
        r2 = outdir / "reads_trimmed/{sample}/{sample}_R2.trimmed.fastq.gz"
    output:
        outdir / "reads_trimmed/{sample}/{sample}_R1.trimmed_fastqc.html",
        outdir / "reads_trimmed/{sample}/{sample}_R2.trimmed_fastqc.html"
    log:
        temp(outdir / "reads_trimmed/{sample}.fastqc.log")
    threads:
        2
    shell:
        """
            fastqc \
                -t {threads} \
                -o "{outdir}/reads_trimmed/{wildcards.sample}" \
                {input.r1} {input.r2} &> {log}
        """

rule fastqc_pretrim_se:
    input:
        r1 = outdir / "reads_raw/{sample}/{sample}_R1.fastq.gz",
    output:
        html = outdir / "reads_raw/{sample}/{sample}_R1_fastqc.html"
    log:
        temp(outdir / "reads_raw/{sample}.fastqc.log")
    threads:
        1
    shell:
        """
            fastqc \
                -t {threads} \
                -o "{outdir}/reads_raw/{wildcards.sample}" \
                {input.r1} &> {log}
        """


rule fastqc_posttrim_se:
    input:
        r1 = outdir / "reads_trimmed/{sample}/{sample}_R1.trimmed.fastq.gz",
    output:
        html = outdir / "reads_trimmed/{sample}/{sample}_R1.trimmed_fastqc.html"
    log:
        temp(outdir / "reads_trimmed/{sample}.fastqc.log")
    threads:
        1
    shell:
        """
            fastqc \
                -t {threads} \
                -o "{outdir}/reads_trimmed/{wildcards.sample}" \
                {input.r1} &> {log}
        """

        


#            FastQC - A high throughput sequence QC analysis tool
#
# SYNOPSIS
# 
# 	fastqc seqfile1 seqfile2 .. seqfileN
# 
#     fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
#            [-c contaminant file] seqfile1 .. seqfileN
# 
# DESCRIPTION
# 
#     FastQC reads a set of sequence files and produces from each one a quality
#     control report consisting of a number of different modules, each one of 
#     which will help to identify a different potential type of problem in your
#     data.
#     
#     If no files to process are specified on the command line then the program
#     will start as an interactive graphical application.  If files are provided
#     on the command line then the program will run with no user interaction
#     required.  In this mode it is suitable for inclusion into a standardised
#     analysis pipeline.
#     
#     The options for the program as as follows:
#     
#     -h --help       Print this help file and exit
#     
#     -v --version    Print the version of the program and exit
#     
#     -o --outdir     Create all output files in the specified output directory.
#                     Please note that this directory must exist as the program
#                     will not create it.  If this option is not set then the 
#                     output file for each sequence file is created in the same
#                     directory as the sequence file which was processed.
#                     
#     --casava        Files come from raw casava output. Files in the same sample
#                     group (differing only by the group number) will be analysed
#                     as a set rather than individually. Sequences with the filter
#                     flag set in the header will be excluded from the analysis.
#                     Files must have the same names given to them by casava
#                     (including being gzipped and ending with .gz) otherwise they
#                     won't be grouped together correctly.
#                     
#     --nano          Files come from nanopore sequences and are in fast5 format. In
#                     this mode you can pass in directories to process and the program
#                     will take in all fast5 files within those directories and produce
#                     a single output file from the sequences found in all files.                    
#                     
#     --nofilter      If running with --casava then don't remove read flagged by
#                     casava as poor quality when performing the QC analysis.
#                    
#     --extract       If set then the zipped output file will be uncompressed in
#                     the same directory after it has been created.  By default
#                     this option will be set if fastqc is run in non-interactive
#                     mode.
#                     
#     -j --java       Provides the full path to the java binary you want to use to
#                     launch fastqc. If not supplied then java is assumed to be in
#                     your path.
#                    
#     --noextract     Do not uncompress the output file after creating it.  You
#                     should set this option if you do not wish to uncompress
#                     the output when running in non-interactive mode.
#                     
#     --nogroup       Disable grouping of bases for reads >50bp. All reports will
#                     show data for every base in the read.  WARNING: Using this
#                     option will cause fastqc to crash and burn if you use it on
#                     really long reads, and your plots may end up a ridiculous size.
#                     You have been warned!
#                     
#     --min_length    Sets an artificial lower limit on the length of the sequence
#                     to be shown in the report.  As long as you set this to a value
#                     greater or equal to your longest read length then this will be
#                     the sequence length used to create your read groups.  This can
#                     be useful for making directly comaparable statistics from 
#                     datasets with somewhat variable read lengths.
#                     
#     -f --format     Bypasses the normal sequence file format detection and
#                     forces the program to use the specified format.  Valid
#                     formats are bam,sam,bam_mapped,sam_mapped and fastq
#                     
#     -t --threads    Specifies the number of files which can be processed
#                     simultaneously.  Each thread will be allocated 250MB of
#                     memory so you shouldn't run more threads than your
#                     available memory will cope with, and not more than
#                     6 threads on a 32 bit machine
#                   
#     -c              Specifies a non-default file which contains the list of
#     --contaminants  contaminants to screen overrepresented sequences against.
#                     The file must contain sets of named contaminants in the
#                     form name[tab]sequence.  Lines prefixed with a hash will
#                     be ignored.
# 
#     -a              Specifies a non-default file which contains the list of
#     --adapters      adapter sequences which will be explicity searched against
#                     the library. The file must contain sets of named adapters
#                     in the form name[tab]sequence.  Lines prefixed with a hash
#                     will be ignored.
#                     
#     -l              Specifies a non-default file which contains a set of criteria
#     --limits        which will be used to determine the warn/error limits for the
#                     various modules.  This file can also be used to selectively 
#                     remove some modules from the output all together.  The format
#                     needs to mirror the default limits.txt file found in the
#                     Configuration folder.
#                     
#    -k --kmers       Specifies the length of Kmer to look for in the Kmer content
#                     module. Specified Kmer length must be between 2 and 10. Default
#                     length is 7 if not specified.
#                     
#    -q --quiet       Supress all progress messages on stdout and only report errors.
#    
#    -d --dir         Selects a directory to be used for temporary files written when
#                     generating report images. Defaults to system temp directory if
#                     not specified.
#                     
# BUGS
# 
#     Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
#     or in www.bioinformatics.babraham.ac.uk/bugzilla/
#                    
#     
