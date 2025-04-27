#!/usr/bin/env python3

import os
import argparse
import subprocess
from functools import partial
import sys

print_warning = partial(print, "WARNING: ", sep="", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description='Wrapper for featureCounts tool.')

    # Mandatory arguments
    parser.add_argument('-a', '--annotation', required=True, type=str, help="Name of the annotation file.")
    parser.add_argument('-o', '--output', required=True, type=str, help="Name of output file including read counts.")
    parser.add_argument('-i', '--input_files', nargs='+', help="A list of SAM or BAM format files.", required=True)

    # Strandness
    parser.add_argument('-s', '--strand', type=str, help="Perform strand-specific read counting.", required=False)


    # Optional arguments
    parser.add_argument('-F', '--annotation_format', type=str, default='GTF', help="Format of the provided annotation file (GTF or SAF).", choices=['GTF', 'SAF', "BED"])
    parser.add_argument('-t', '--feature_type', type=str, help="Feature type(s) in a GTF annotation.", required=False)
    parser.add_argument('-g', '--attribute_type', type=str, help="Attribute type in GTF annotation.", required=False)
    parser.add_argument('--extraAttributes', type=str, help="Extract extra attribute types from GTF annotation.")
    parser.add_argument('-A', '--alias_file', type=str, help="Provide a chromosome name alias file.")
    parser.add_argument("--names", nargs='+', help="Provide a list of sample names for the input files.")
    parser.add_argument("--bed-key", choices=("locus", "name"), default="locus", )
    parser.add_argument("--gene-table", type=str, default=None, help="Gene table file for gene name conversion.")

    # Level of summarization
    parser.add_argument('-f', '--feature_level', action='store_true', help="Perform read counting at feature level.")

    # Overlap between reads and features
    parser.add_argument('-O', '--overlap', action='store_true', help="Assign reads to all overlapping meta-features.")
    parser.add_argument('--minOverlap', type=int, help="Minimum number of overlapping bases in a read required for assignment.")
    parser.add_argument('--fracOverlap', type=float, help="Minimum fraction of overlapping bases in a read required for assignment.")
    parser.add_argument('--fracOverlapFeature', type=float, help="Minimum fraction of overlapping bases in a feature required for assignment.")
    parser.add_argument('--largestOverlap', action='store_true', help="Assign reads to the feature with the largest number of overlapping bases.")
    parser.add_argument('--nonOverlap', type=int, help="Maximum number of non-overlapping bases in a read allowed when assigning to a feature.")
    parser.add_argument('--nonOverlapFeature', type=int, help="Maximum number of non-overlapping bases in a feature allowed in read assignment.")
    parser.add_argument('--readExtension5', type=int, help="Extend reads upstream by specified bases from 5' end.")
    parser.add_argument('--readExtension3', type=int, help="Extend reads upstream by specified bases from 3' end.")
    parser.add_argument('--read2pos', type=str, choices=['5', '3'], help="Reduce reads to their 5' or 3' most base.")

    # Multi-mapping reads
    parser.add_argument('-M', '--multi_mapping', action='store_true', help="Count multi-mapping reads.")

    # Fractional counting
    parser.add_argument('--fraction', action='store_true', help="Assign fractional counts to features (use with -M or -O).")

    # Read filtering
    parser.add_argument('-Q', '--min_quality', type=int, help="Minimum mapping quality score a read must satisfy to be counted.")
    parser.add_argument('--splitOnly', action='store_true', help="Count split alignments only.")
    parser.add_argument('--nonSplitOnly', action='store_true', help="Count non-split alignments only.")
    parser.add_argument('--primary', action='store_true', help="Count primary alignments only.")
    parser.add_argument('--ignoreDup', action='store_true', help="Ignore duplicate reads in read counting.")

    # Exon-exon junctions
    parser.add_argument('-J', '--junctions', action='store_true', help="Count number of reads supporting each exon-exon junction.")
    parser.add_argument('-G', '--reference', type=str, help="Provide the reference sequences used in read mapping.")

    # Parameters specific to paired-end reads
    parser.add_argument('-p', '--paired_end', action='store_true', help="Input data contains paired-end reads.")
    parser.add_argument('--countReadPairs', action='store_true', help="Count read pairs instead of reads (for paired-end reads).")
    parser.add_argument('-B', '--both_ends', action='store_true', help="Only count read pairs that have both ends aligned.")
    parser.add_argument('-P', '--valid_distance', action='store_true', help="Check validity of paired-end distance when counting read pairs.")
    parser.add_argument('-d', '--min_length', type=int, help="Minimum fragment length (default: 50).")
    parser.add_argument('-D', '--max_length', type=int, help="Maximum fragment length (default: 600).")
    parser.add_argument('-C', '--same_chromosome', action='store_true', help="Do not count read pairs with ends mapping to different chromosomes.")
    parser.add_argument('--donotsort', action='store_true', help="Do not sort reads in BAM/SAM input.")

    # Number of CPU threads
    parser.add_argument('-T', '--threads', type=int, default=min(8, os.cpu_count()), help="Number of threads (default: {}).".format(min(8, os.cpu_count())))

    # Read groups
    parser.add_argument('--byReadGroup', action='store_true', help="Assign reads by read group (RG tag required).")

    # Long reads
    parser.add_argument('-L', '--long_reads', action='store_true', help="Count long reads such as Nanopore and PacBio reads.")

    # Assignment results for each read
    parser.add_argument('-R', '--detailed_results', type=str, choices=['CORE', 'SAM', 'BAM'], help="Output detailed assignment results for each read.")
    parser.add_argument('--Rpath', type=str, help="Specify directory to save detailed assignment results.")

    # Miscellaneous
    parser.add_argument('--tmpDir', type=str, help="Directory to save intermediate files.")
    parser.add_argument('--maxMOp', type=int, help="Maximum number of M operations allowed in CIGAR string (default: 10).")
    parser.add_argument('--verbose', action='store_true', help="Output verbose information for debugging.")
    parser.add_argument('-v', '--version', action='store_true', help="Output version of the program.")

    args = parser.parse_args()


    assert args.strand is not None, "Strandness must be provided (-s 0/1/2)."

    if args.annotation_format == "GTF":
        assert args.feature_type, "Feature type must be provided for GTF annotation."
        assert args.attribute_type, "Attribute type must be provided for GTF annotation."
    elif args.annotation_format == "BED":
        tmp_bed = args.annotation + '.tmp.saf'
        if args.bed_key == "locus":
            os.system("bed2saf.py {} > {}".format(args.annotation, tmp_bed))
        else:
            os.system("bed2saf.py -n {} > {}".format(args.annotation, tmp_bed))
        args.annotation = tmp_bed
        args.annotation_format = "SAF"

    # Build the featureCounts command
    cmd = ['featureCounts']

    if args.names:
        assert len(args.names) == len(args.input_files), "Number of sample names must match the number of input files."
        assert len(set(args.names)) == len(args.names), "Sample names must be unique."

    if args.annotation == "hg38":
        hg38_default = "/home/kenchen/db/gencode/GRCh38/v46/gencode.v46.annotation.gtf"
        print_warning("Using default annotation file for human genome (GRCh38): " + hg38_default)
        args.annotation = hg38_default


    # Mandatory arguments
    cmd += ['-a', args.annotation, '-o', args.output]
    cmd += args.input_files


    # Optional arguments
    if args.annotation_format:
        cmd += ['-F', args.annotation_format]
    if args.feature_type:
        cmd += ['-t', args.feature_type]
    if args.attribute_type:
        cmd += ['-g', args.attribute_type]
    if args.extraAttributes:
        cmd += ['--extraAttributes', args.extraAttributes]
    if args.alias_file:
        cmd += ['-A', args.alias_file]
    if args.feature_level:
        cmd.append('-f')
    if args.overlap:
        cmd.append('-O')
    if args.minOverlap:
        cmd += ['--minOverlap', str(args.minOverlap)]
    if args.fracOverlap:
        cmd += ['--fracOverlap', str(args.fracOverlap)]
    if args.fracOverlapFeature:
        cmd += ['--fracOverlapFeature', str(args.fracOverlapFeature)]
    if args.largestOverlap:
        cmd.append('--largestOverlap')
    if args.nonOverlap:
        cmd += ['--nonOverlap', str(args.nonOverlap)]
    if args.nonOverlapFeature:
        cmd += ['--nonOverlapFeature', str(args.nonOverlapFeature)]
    if args.readExtension5:
        cmd += ['--readExtension5', str(args.readExtension5)]
    if args.readExtension3:
        cmd += ['--readExtension3', str(args.readExtension3)]
    if args.read2pos:
        cmd += ['--read2pos', args.read2pos]
    if args.multi_mapping:
        cmd.append('-M')
    if args.fraction:
        cmd.append('--fraction')
    if args.min_quality:
        cmd += ['-Q', str(args.min_quality)]
    if args.splitOnly:
        cmd.append('--splitOnly')
    if args.nonSplitOnly:
        cmd.append('--nonSplitOnly')
    if args.primary:
        cmd.append('--primary')
    if args.ignoreDup:
        cmd.append('--ignoreDup')
    if args.strand:
        cmd += ['-s', args.strand]
    if args.junctions:
        cmd.append('-J')
    if args.reference:
        cmd += ['-G', args.reference]
    if args.paired_end:
        cmd.append('-p')
    if args.countReadPairs:
        cmd.append('--countReadPairs')
    if args.both_ends:
        cmd.append('-B')
    if args.valid_distance:
        cmd.append('-P')
    if args.min_length:
        cmd += ['-d', str(args.min_length)]
    if args.max_length:
        cmd += ['-D', str(args.max_length)]
    if args.same_chromosome:
        cmd.append('-C')
    if args.donotsort:
        cmd.append('--donotsort')
    if args.threads:
        cmd += ['-T', str(args.threads)]
    if args.byReadGroup:
        cmd.append('--byReadGroup')
    if args.long_reads:
        cmd.append('-L')
    if args.detailed_results:
        cmd += ['-R', args.detailed_results]
    if args.Rpath:
        cmd += ['--Rpath', args.Rpath]
    if args.tmpDir:
        cmd += ['--tmpDir', args.tmpDir]
    if args.maxMOp:
        cmd += ['--maxMOp', str(args.maxMOp)]
    if args.verbose:
        cmd.append('--verbose')
    if args.version:
        cmd.append('-v')
    
    print("#Command: ", ' '.join(cmd), file=sys.stderr)
    # Run the command
    subprocess.run(cmd)

    # rename sample names
    if args.names:
        with open(args.output, 'r') as f, open(args.output + '.tmp', 'w') as g:
            g.write("##convert raw names to: {}".format(' '.join(args.names)) + '\n')
            for line in f:
                if line.startswith('#'):
                    g.write(line)
                elif line.startswith('Geneid'):
                    fields = line.rstrip('\n').split('\t')
                    fields = fields[:-len(args.input_files)] + args.names
                    g.write('\t'.join(fields) + '\n')
                else:
                    g.write(line)
        os.system("mv {} {}".format(args.output + '.tmp', args.output))
    if args.gene_table:
        gt = pd.read_table(args.gene_table, index_col=0)


if __name__ == '__main__':
    main()

