"""Compare PCA compressibility: raw expression data (68 tissues x 8294 genes with
TPM>=1 in every tissue, log10, centered-only) vs each model's fitted parameter space
(68 tissues x theta, standardized). Both are single global PCAs; the resulting top-3
%variance-explained numbers are global (one triple per PCA), tabulated identically on
every tissue row so the whole comparison lives in one tidy, joinable CSV. Per model,
counts how many of its top-3 PCs (rank-matched: PC1 vs PC1, PC2 vs PC2, PC3 vs PC3)
explain more variance than the expression data's top-3 PCs."""
import sys
import pandas as pd
import numpy as np
from mixture_core import theta_from_fit, canonicalize_by_weight, unflatten_fit
from model_registry import MODELS

IN_EXPR_MATRIX = sys.argv[1]
MODEL_FITS_DIR = sys.argv[2]
OUT_CSV = sys.argv[3]


def top3_var_ratio(X, standardize):
    if standardize:
        mean = X.mean(axis=0)
        std = X.std(axis=0, ddof=1)
        std[std < 1e-8] = 1.0
        Xz = (X - mean) / std
    else:
        Xz = X - X.mean(axis=0)
    _, S, _ = np.linalg.svd(Xz, full_matrices=False)
    var_ratio = (S ** 2) / (S ** 2).sum()
    return var_ratio[:3]


expr = pd.read_parquet(IN_EXPR_MATRIX)
tissues = [c for c in expr.columns if c != "gene"]
X_expr = expr[tissues].values.T  # tissues x genes
expr_var = top3_var_ratio(X_expr, standardize=False)

model_vars = {}
for model, comp_types in MODELS.items():
    fits = pd.read_csv(f"{MODEL_FITS_DIR}/{model}_fits.csv")
    X_model = np.stack([
        theta_from_fit(comp_types, *canonicalize_by_weight(
            comp_types, *unflatten_fit(comp_types, row)))
        for _, row in fits.iterrows()
    ])
    model_vars[model] = top3_var_ratio(X_model, standardize=True)

rows = []
for t in tissues:
    row = {"tissue": t, "expr_pc1_var": expr_var[0], "expr_pc2_var": expr_var[1], "expr_pc3_var": expr_var[2]}
    for model, mvar in model_vars.items():
        row[f"{model}_beat_count"] = int(sum(mvar[i] > expr_var[i] for i in range(3)))
    rows.append(row)

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)

print("expr top3 var:", expr_var, file=sys.stderr)
for model, mvar in model_vars.items():
    count = int(sum(mvar[i] > expr_var[i] for i in range(3)))
    print(f"{model:22s} top3 var: {mvar} beat_count: {count}", file=sys.stderr)
