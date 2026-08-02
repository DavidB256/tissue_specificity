"""Combine per-model fit tables into one wide comparison CSV: rows=tissue, cols=aic/bic/ks per model."""
import sys
import pandas as pd
from model_registry import MODELS

RESULTS_DIR = sys.argv[1]
OUT_CSV = sys.argv[2]

wide = None
for model in MODELS:
    df = pd.read_csv(f"{RESULTS_DIR}/{model}_fits.csv")[["tissue", "aic", "bic", "ks_stat"]]
    df = df.rename(columns={"aic": f"{model}_aic", "bic": f"{model}_bic", "ks_stat": f"{model}_ks"})
    wide = df if wide is None else wide.merge(df, on="tissue")
wide.to_csv(OUT_CSV, index=False)
