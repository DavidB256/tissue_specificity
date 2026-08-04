"""One-time cleanup: clip comp{i}_a columns (skewnorm shape param) to [-100, 100] in
existing fit CSVs. Not part of the Snakefile pipeline -- run manually, once, against
each model_fits/ dir. Mirrors the a-clip now enforced during fitting in mixture_core.py."""
import sys
import glob
import pandas as pd
from model_registry import MODELS

DATA_DIRS = sys.argv[1:]


def skew_a_cols(comp_types):
    return [f"comp{i + 1}_a" for i, ct in enumerate(comp_types) if ct == "skewnorm"]


def clip_file(path, a_cols):
    df = pd.read_csv(path)
    cols = [c for c in a_cols if c in df.columns]
    if not cols:
        return False
    changed = False
    for c in cols:
        n = ((df[c] > 100) | (df[c] < -100)).sum()
        if n:
            changed = True
            df[c] = df[c].clip(-100, 100)
    if changed:
        df.to_csv(path, index=False)
    return changed


for d in DATA_DIRS:
    for model, comp_types in MODELS.items():
        a_cols = skew_a_cols(comp_types)
        if not a_cols:
            continue
        for path in sorted(glob.glob(f"{d}/{model}_fits.csv")):
            changed = clip_file(path, a_cols)
            print(f"{path}: {'clipped' if changed else 'no change'}", file=sys.stderr)
        for subdir in ("by_model_tissue", "by_condition"):
            per_item_files = sorted(glob.glob(f"{d}/{subdir}/{model}/*.csv"))
            n_clipped = sum(clip_file(p, a_cols) for p in per_item_files)
            if per_item_files:
                print(f"{d}/{subdir}/{model}/: clipped {n_clipped}/{len(per_item_files)} files", file=sys.stderr)
