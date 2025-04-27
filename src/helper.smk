#!/usr/bin/env snakemake

## Author: Ken Chen <chenkenbio@gmail.com>
## Date: 2025-03-18
import pandas as pd
from typing import List, Dict, Any, Tuple
import warnings

def load_sample_info(fn: str, allow_multi_runs: bool=False) -> Dict[str, Any]:
    sample_dict = dict()
    kv = dict()
    with open(fn, 'r') as infile:
        for l in infile:
            if l.startswith("##"):
                l = l[2:]
                assert '=' in l
                key, value = l.strip().split('=')
                key, value = key.strip(), value.strip()
                kv[key] = value
    assert "path" in kv
    df = pd.read_table(fn, comment="#")
    for i, row in df.iterrows():
        row = row.to_dict()
        sample = row["sample"]
        row["R1"] = os.path.join(kv["path"], row["R1"])
        if "R2" in row:
            row["R2"] = os.path.join(kv["path"], row["R2"])
        if "genome" not in row and "genome" in kv:
            row["genome"] = kv["genome"]

        if "library_type" not in row and "library_type" in kv:
            assert kv["library_type"] in ["RF", "FR", "R", "F", "UNSTRANDED"]
            row["library_type"] = kv["library_type"]

        if sample in sample_dict:
            if allow_multi_runs:
                warnings.warn(f"Duplicate sample name {sample} in {fn}.")
                for k, v1 in row.items():
                    v0 = sample_dict[sample][k]
                    if k in ["R1", "R2"]:
                        if isinstance(v0, str):
                            row[k] = [v0] + [v1]
                        else:
                            assert isinstance(v0, list), f"Sample {sample} has different {k} in {fn}: {v0} vs {v1}"
                            row[k] = v0 + [v1]
                    else:
                        assert v0 == v1, f"Sample {sample} has different {k} in {fn}: {v0} vs {v1}"
            else:
                raise ValueError(f"Duplicate sample name {sample} in {fn}.")

        sample_dict[sample] = row
            
    return sample_dict



