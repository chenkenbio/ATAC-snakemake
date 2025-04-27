#!/usr/bin/env python3

import warnings
import os
import sys
import pandas as pd
import argparse
import numpy as np

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--featureCounts', required=True, help='featureCounts count file', nargs='+')
    p.add_argument("-l", "--label", nargs='+', help="labels for each featureCounts file")
    p.add_argument("-p", "--prefix", required=True, help="prefix for output file")
    p.add_argument("--ignore-summary", action='store_true', help="Ignore summary file")
    p.add_argument("--no-header", action='store_true', help="ignore header in featureCounts file")
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--annotation", "-a", help="Annotation file for feature id")
    return p

if __name__ == '__main__':
    args = parse_args().parse_args()

    for fn in args.featureCounts:
        assert os.path.exists(fn), f"File not found: {fn}"
        if not args.ignore_summary:
            assert os.path.exists(fn + '.summary'), f"File not found: {fn}.summary"
    if args.label is not None:
        assert len(args.label) == len(args.featureCounts), f"Number of labels does not match number of featureCounts files: {len(args.label)} != {len(args.featureCounts)}"
        for name, fn in zip(args.label, args.featureCounts):
            if name not in fn:
                warnings.warn(f"Warning: label {name} not found in filename {fn}")

    counts = list()
    summary = list()
    index = None
    headers = list()
    for i, fn in enumerate(args.featureCounts):
        print(f"Reading {fn}")
        if args.verbose:
            print(pd.read_table(fn, delimiter='\t', index_col=0, comment='#').head())
        with open(fn) as f:
            header = f.readline().rstrip('\n')
            while header.startswith("##"):
                header = f.readline().rstrip('\n')
            assert header.startswith("# Program:featureCounts"), f"File {fn} is not a featureCounts file, first line: {header}"
            headers.append(header)
        df = pd.read_table(fn, delimiter='\t', index_col=0, comment='#')
        if args.label is not None:
            cols = list(df.columns)
            cols[-1] = args.label[i]
            df.columns = cols
        if index is None:
            index = df.index
            annotation = df.iloc[:, range(0, 5)]
        else:
            assert sum([a == b for a, b in zip(index, df.index)]) == len(index), f"Index mismatch: {index} != {df.index}"
        counts.append(df.iloc[:, 5:])
        if not args.ignore_summary:
            df = pd.read_table(fn + '.summary', delimiter='\t', index_col=0, comment='#')
            if args.label is not None:
                cols = list(df.columns)
                cols[-1] = args.label[i]
                df.columns = cols
            summary.append(df)

    # write header before writing counts
    counts = pd.concat([annotation] + counts, axis=1)
    summary = pd.concat(summary, axis=1)
    if args.prefix.endswith('.counts'):
        args.prefix = args.prefix[:-7]
    if not args.no_header:
        with open(args.prefix + '.counts', 'w') as f:
            f.write('\n'.join(headers) + '\n')
        counts.to_csv(args.prefix + '.counts', sep='\t', mode='a')
    else:
        counts.to_csv(args.prefix + '.counts', sep='\t')
    summary.to_csv(args.prefix + '.counts.summary', sep='\t')




