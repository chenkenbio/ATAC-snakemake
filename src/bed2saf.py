#!/usr/bin/env python3
"""
Author: Ken Chen (https://github.com/chenkenbio)
Date: 2024-09-25
"""

import argparse
import os
import sys
import math
import warnings
import numpy as np
from biock2 import auto_open

def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("bed")
    p.add_argument("-n", "--name", action='store_true', help="Use the 4th column as name")
    p.add_argument("--add-prefix", help="Add prefix to the name", default=None)
    p.add_argument("--allow-dup", action='store_true', help="Allow duplicate names")
    # p.add_argument('--seed', type=int, default=2020)
    return p


if __name__ == "__main__":
    args = get_args().parse_args()

    seen = set()

    with auto_open(args.bed, 'rt') as infile:
        cnt = 0
        for l in infile:
            if l.startswith("#"):
                continue
            fields = l.strip().split("\t")
            chrom, start, end = fields[:3]
            if len(fields) >= 6:
                strand = fields[5]
            else:
                strand = "."
            if args.name:
                name = fields[3]
            else:
                if args.add_prefix:
                    name = f"{args.add_prefix}{cnt}_{chrom}_{start}_{end}_{strand}"
                else:
                    name = f"{chrom}_{start}_{end}_{strand}"
                    if name in seen:
                        continue
            cnt += 1
            if not args.allow_dup:
                assert name not in seen, f"{name} already exists"
            seen.add(name)
            # print(f"{chrom}\t{start}\t{end}\t{name}\t.\t{strand}")
            print(f"{name}\t{chrom}\t{int(start)+1}\t{end}\t{strand}")
