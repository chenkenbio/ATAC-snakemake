#!/usr/bin/env python3
"""
Author: Ken Chen (https://github.com/chenkenbio)
"""

import argparse
import os
import sys
import gzip
import math
import json
import time
import string
import random
import warnings
import numpy as np
import pandas as pd
from tqdm import tqdm
from glob import glob
from collections import defaultdict, OrderedDict
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Union, Callable
from subprocess import Popen, PIPE

def auto_open(filename, mode='r', **kwargs):
    """
    Open a file automatically based on its extension.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode=mode, **kwargs)
    elif filename.endswith('.bz2'):
        import bz2
        return bz2.BZ2File(filename, mode=mode, **kwargs)
    elif filename.endswith('.xz'):
        import lzma
        return lzma.open(filename, mode=mode, **kwargs)
    else:
        return open(filename, mode=mode, **kwargs)


def random_string(n):
    random.seed(time.time() % 3533)
    return ''.join(random.choices(string.ascii_letters + string.digits, k=n))


def random_prefix(basename="tmpfile", tmp_dir: str=None) -> str:
    if tmp_dir is None:
        tmp_dir = os.environ.get("TMPDIR", os.path.expanduser("~/tmp"))
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    prefix = os.path.join(
        tmp_dir, 
        "{}_{}".format(basename, random_string(8))
    )
    while len(glob(prefix + "*")) > 0:
        prefix = os.path.join(
            tmp_dir, 
            "{}_{}".format(basename, random_string(8))
        )
    return prefix


def run_bash(cmd, check_rc: bool=True) -> Tuple[int, str, str]:
    r"""
    Return
    -------
    rc : return code
    out : output
    err : error
    """
    p = Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    out, err = out.decode('utf8'), err.decode('utf8')
    rc = p.returncode
    if check_rc and rc != 0:
        raise RuntimeError("Command {} failed with return code {}\n{}".format(cmd, rc, err))
    return (rc, out, err)


def strip_end(s: str, suffix: str) -> str:
    if s.endswith(suffix):
        return s[:-len(suffix)]
    else:
        return s



def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--bam", '-i', help="input bam file")
    p.add_argument("-o", "--output", help="output bam file", required=False)
    p.add_argument("--keep-unmapped", action="store_true", help="keep unmapped reads")
    p.add_argument("--rmdup", action="store_true", help="remove duplicates")
    p.add_argument("--min-mapq", "-mq", type=int, default=0, help="minimum mapping quality")
    p.add_argument("-t", "--threads", type=int, default=8, help="number of threads")
    p.add_argument("--blacklist", "-b", help="blacklist file", default=None)
    p.add_argument("--exclude-chromosomes", "-e", help="exclude chromosomes", nargs='+', default=None)
    p.add_argument("--chr-only", action="store_true", help="only keep chromosomes started with chr")
    # p.add_argument('--seed', type=int, default=2020)
    return p

def get_chromomes_from_bam_header(bam_file) -> Dict[str, int]:
    rc, out, err = run_bash(f"samtools view -H {bam_file} | grep '^@SQ'")
    assert rc == 0, f"Error getting chromosomes from bam header: {err}"
    chromsize = dict()
    for line in out.split('\n'):
        if line.startswith('@SQ'):
            parts = line.split('\t')
            chrom = parts[1].split(':')[1]
            size = int(parts[2].split(':')[1])
            chromsize[chrom] = size
    return chromsize

def filtering_bam(
        bam, 
        output, 
        *, 
        min_mapq: int=0, 
        keep_unmapped: bool,
        rmdup: bool, 
        chr_only: bool,
        blacklist_bed, 
        threads: int=8,
        exclude_chromosomes: List[str],
    ) -> str:

    if blacklist_bed is not None or exclude_chromosomes is not None or chr_only:
        chromsizes = get_chromomes_from_bam_header(bam)

        selected_bed = random_prefix() + ".bed"
        with open(selected_bed, 'w') as outfile:
            for chrom, size in chromsizes.items():
                if chr_only and not chrom.startswith("chr"):
                    continue
                outfile.write(f"{chrom}\t0\t{size}\n")

        if exclude_chromosomes or blacklist_bed:
            exclude_bed = random_prefix() + ".bed"
            with open(exclude_bed, 'w') as outfile:
                if blacklist_bed is not None:
                    with auto_open(blacklist_bed, 'rt') as infile:
                        for l in infile:
                            chrom, start, end = l.rstrip().split('\t')[:3]
                            outfile.write(f"{chrom}\t{start}\t{end}\n")
                if exclude_chromosomes:
                    for c in exclude_chromosomes:
                        if c in chromsizes:
                            size = chromsizes[c]
                            outfile.write(f"{c}\t0\t{size}\n")
                        else:
                            print(f"Warning: chromosome {c} not found in bam header")
            os.system(f"bedtools subtract -a {selected_bed} -b {exclude_bed} > {selected_bed}.tmp && mv {selected_bed}.tmp {selected_bed} && rm {exclude_bed}")
    else:
        selected_bed = None

    cmd = ' '.join([
        f"samtools view -b {bam}", 
        f"-@ {threads}",
        f"-F 4" if not keep_unmapped else "",
        f"-F 1024" if rmdup else "",
        f"-q {min_mapq}" if min_mapq > 0 else "",
        f"-L {selected_bed}" if selected_bed else "",
        f"> {output}"
    ])
    print(f"Running command: {cmd}", file=sys.stderr)
    os.system(cmd)
    if selected_bed:
        os.remove(selected_bed)


if __name__ == "__main__":
    args = get_args().parse_args()
    
    output = args.output if args.output else (strip_end(args.bam, ".bam") + ".filtered.bam")
    filtering_bam(
        bam=args.bam,
        output=output,
        min_mapq=args.min_mapq,
        keep_unmapped=args.keep_unmapped,
        rmdup=args.rmdup,
        chr_only=args.chr_only,
        blacklist_bed=args.blacklist,
        threads=args.threads,
        exclude_chromosomes=args.exclude_chromosomes,
    )
