"""One-time migration: rewrite existing {model}_fits.csv files, and their
by_model_tissue/{model}/{tissue}.csv sources, from JSON-encoded weights/params
columns to flat one-value-per-cell columns (weight_1.., comp1_<param>..). Not part
of the Snakefile pipeline -- run manually, once, against each model_fits/ dir."""
import sys
import os
import json
import glob
import pandas as pd
from mixture_core import flatten_fit
from model_registry import MODELS

DATA_DIRS = sys.argv[1:]


def migrate_file(path, comp_types):
    df = pd.read_csv(path)
    if "weights" not in df.columns or "params" not in df.columns:
        return False
    rows = []
    for _, row in df.iterrows():
        weights = json.loads(row["weights"])
        params = json.loads(row["params"])
        flat = flatten_fit(comp_types, weights, params)
        base = row.drop(labels=["weights", "params"]).to_dict()
        base.update(flat)
        rows.append(base)
    pd.DataFrame(rows).to_csv(path, index=False)
    return True


for d in DATA_DIRS:
    for path in sorted(glob.glob(f"{d}/*_fits.csv")):
        model = os.path.basename(path)[: -len("_fits.csv")]
        if model not in MODELS:
            print(f"{path}: skipped (not a registered model, leaving as-is)", file=sys.stderr)
            continue
        changed = migrate_file(path, MODELS[model])
        print(f"{path}: {'migrated' if changed else 'skipped (already flat)'}", file=sys.stderr)

    for subdir in ("by_model_tissue", "by_condition"):
        for model, comp_types in MODELS.items():
            per_item_files = sorted(glob.glob(f"{d}/{subdir}/{model}/*.csv"))
            n_migrated = sum(migrate_file(p, comp_types) for p in per_item_files)
            if per_item_files:
                print(f"{d}/{subdir}/{model}/: migrated {n_migrated}/{len(per_item_files)} files", file=sys.stderr)
