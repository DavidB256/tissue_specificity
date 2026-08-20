"""plot_tau_range_unweighted_vs_ridge.py restricted to "liver" genes: max-expressing
tissue (all-tissue, TPM floored) is a Liver_* tissue (same logic as
gls_tau.gene_subset_mask's "maxneural", but for Liver_* instead of Brain_*)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_utils import annotate_diag_stats
from gls_tau import floor_tpm

IN_GCT = sys.argv[1]
IN_UNWEIGHTED = sys.argv[2]
IN_RIDGE = sys.argv[3]
IN_TAU = sys.argv[4]
OUT_PNG = sys.argv[5]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = floor_tpm(gct[tissues].values.astype(float))
liver_idx = np.array([i for i, t in enumerate(tissues) if "Liver" in t])
is_liver = np.isin(mat.argmax(axis=1), liver_idx)
liver = pd.DataFrame({"gene": gct["Name"].values, "is_liver": is_liver})

u = pd.read_csv(IN_UNWEIGHTED)[["gene", "tau_range"]].rename(columns={"tau_range": "tau_range_unweighted"})
r = pd.read_csv(IN_RIDGE)[["gene", "tau_range"]].rename(columns={"tau_range": "tau_range_ridge"})
tau = pd.read_csv(IN_TAU)[["gene", "tau_unweighted"]]
df = u.merge(r, on="gene").merge(tau, on="gene").merge(liver, on="gene").dropna()
df = df[df["is_liver"]]

fig, ax = plt.subplots(figsize=(6, 6))
sc = ax.scatter(df["tau_range_unweighted"], df["tau_range_ridge"], c=df["tau_unweighted"],
                 cmap="viridis", s=4, alpha=0.1)
fig.colorbar(sc, ax=ax, label="tau_unweighted (all-tissue)")
hi = max(df["tau_range_unweighted"].max(), df["tau_range_ridge"].max())
ax.plot([0, hi], [0, hi], "--", color="gray", lw=1)
ax.set_xlabel("tau_range (unweighted, jackknife)")
ax.set_ylabel("tau_range (gls_ridge, jackknife)")
ax.set_title(f"jackknife tau_range: unweighted vs ridge, liver genes, n={len(df)} genes")
annotate_diag_stats(ax, df["tau_range_unweighted"], df["tau_range_ridge"])
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
