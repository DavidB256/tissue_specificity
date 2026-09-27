"""Top-20-per-scenario table of genes with the lowest DeltaDeltaTau
(= |delta_tau_w| - |delta_tau|, most negative = free weighting most shrank the tau shift
from collapsing tissue resolution), for onebrain and simplified_0.1. Reports each gene's
top-2 most-expressed tissues (from the full 68-tissue GCT) plus tau/tau_weighted in
the full and alternate (collapsed) scenarios."""
import sys
import numpy as np
import pandas as pd
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_ONEBRAIN_FREE_TAU = sys.argv[2]
IN_ONEBRAIN_REF_TAU = sys.argv[3]
IN_SIMPLIFIED_FREE_TAU = sys.argv[4]
IN_SIMPLIFIED_REF_TAU = sys.argv[5]
OUT_CSV = sys.argv[6]
TOP_N = 20

SCENARIOS = [
    ("onebrain", IN_ONEBRAIN_FREE_TAU, IN_ONEBRAIN_REF_TAU),
    ("simplified_0.1", IN_SIMPLIFIED_FREE_TAU, IN_SIMPLIFIED_REF_TAU),
]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = np.array(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
top2_idx = np.argsort(-mat, axis=1)[:, :2]
top_tissues = pd.DataFrame({
    "gene": gct["Name"].values,
    "gene_symbol": gct["Description"].values,
    "tissue_1": tissues[top2_idx[:, 0]],
    "tissue_2": tissues[top2_idx[:, 1]],
})

tables = []
for scenario, free_tau_csv, ref_tau_csv in SCENARIOS:
    d = pd.read_csv(free_tau_csv)[["gene", "tau_unweighted", "tau_free_weighted"]] \
        .rename(columns={"tau_unweighted": "tau_full", "tau_free_weighted": "tau_weighted_full"})
    r = pd.read_csv(ref_tau_csv)[["gene", "tau_unweighted", "tau_weighted_diag_reg_inv_corr"]] \
        .rename(columns={"tau_unweighted": "tau_alt", "tau_weighted_diag_reg_inv_corr": "tau_weighted_alt"})

    df = d.merge(r, on="gene").merge(top_tissues, on="gene").dropna(
        subset=["tau_full", "tau_weighted_full", "tau_alt"])
    df["delta_tau"] = df["tau_full"] - df["tau_alt"]
    df["delta_tau_w"] = df["tau_weighted_full"] - df["tau_alt"]
    df["delta_delta_tau"] = df["delta_tau_w"].abs() - df["delta_tau"].abs()
    df["scenario"] = scenario

    top = df.nsmallest(TOP_N, "delta_delta_tau")
    tables.append(top[["scenario", "gene", "gene_symbol", "tissue_1", "tissue_2",
                        "tau_full", "tau_weighted_full", "tau_alt", "tau_weighted_alt",
                        "delta_tau", "delta_tau_w", "delta_delta_tau"]])

pd.concat(tables, ignore_index=True).to_csv(OUT_CSV, index=False)
