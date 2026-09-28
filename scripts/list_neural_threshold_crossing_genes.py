"""List "neural" genes (max-expressing tissue, all-tissue floored TPM, is a Brain_*
tissue -- gls_tau.gene_subset_mask's "maxneural") whose tau crosses THRESHOLD going
from tau_unweighted to tau_free_weighted (same crossing convention as
list_crossing_threshold_genes.py). For each crossing gene, reports what fraction of
its total (linear, floor_tpm'd) expression sits in Brain_* tissues, thresholded at
90%/99%/100% -- floor_tpm'd so TPM<1 entries (rounded down to 0 elsewhere in this
pipeline) don't count as expression here either.
Args: IN_GCT IN_TAU_CSV IN_NAME_MAP_CSV THRESHOLD OUT_CSV"""
import sys
import numpy as np
import pandas as pd
from gls_tau import floor_tpm, gene_subset_mask

IN_GCT = sys.argv[1]
IN_TAU_CSV = sys.argv[2]
IN_NAME_MAP_CSV = sys.argv[3]
THRESHOLD = float(sys.argv[4])
OUT_CSV = sys.argv[5]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))

is_neural = gene_subset_mask("maxneural", mat, tissues, None)
brain_idx = np.array([i for i, t in enumerate(tissues) if "Brain" in t])
total = mat.sum(axis=1)
brain_total = mat[:, brain_idx].sum(axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    frac_brain = np.where(total > 0, brain_total / total, np.nan)

genes = pd.DataFrame({
    "gene": gct["Name"].values,
    "is_neural": is_neural,
    "brain_expr_90": frac_brain >= 0.90,
    "brain_expr_99": frac_brain >= 0.99,
    "brain_expr100": frac_brain >= 1.0 - 1e-9,
})

tau = pd.read_csv(IN_TAU_CSV)[["gene", "tau_unweighted", "tau_free_weighted"]].dropna()
df = tau.merge(genes, on="gene")
df = df[df["is_neural"]].copy()

u, w = df["tau_unweighted"], df["tau_free_weighted"]
crossing = ((u > THRESHOLD) & (w < THRESHOLD)) | ((u < THRESHOLD) & (w > THRESHOLD))
df = df[crossing]
df["delta_tau"] = df["tau_free_weighted"] - df["tau_unweighted"]

name_map = pd.read_csv(IN_NAME_MAP_CSV)
name_map["gene_base"] = name_map["gene_id"].str.split(".").str[0]
df["gene_base"] = df["gene"].str.split(".").str[0]
df = df.merge(name_map[["gene_base", "gene_name"]], on="gene_base", how="left")

df = df.sort_values("delta_tau", key=abs, ascending=False)
out = df.rename(columns={"gene": "ensembl_id"})[
    ["ensembl_id", "gene_name", "delta_tau", "brain_expr_90", "brain_expr_99", "brain_expr100"]
]
out.to_csv(OUT_CSV, index=False)

n_unmapped = out["gene_name"].isna().sum()
print(f"neural genes crossing tau={THRESHOLD}: n={len(out)} (unmapped names: {n_unmapped}) -> {OUT_CSV}",
      file=sys.stderr)
