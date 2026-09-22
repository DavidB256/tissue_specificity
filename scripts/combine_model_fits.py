"""Concat per-tissue model-fit CSVs (by_model_tissue/{model}/{tissue}.csv) into one
model_fits/{model}_fits.csv. Args: OUT_CSV, then all per-tissue fit CSVs."""
import sys
import pandas as pd

OUT_CSV = sys.argv[1]
IN_CSVS = sys.argv[2:]

dfs = [pd.read_csv(f) for f in IN_CSVS]
pd.concat(dfs, ignore_index=True).to_csv(OUT_CSV, index=False)
