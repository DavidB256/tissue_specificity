"""Poster fig (free_cluster_weights_max, onebrain leg): histogram of arctanh(tx tau_weighted) - arctanh(gene tau_weighted) over
multi-tx genes (n_tx>=2) of the quantified tx (tx_total>0, i.e. TPM>=1 in >=1 tissue --
the GCT is already floored at TPM<1). Each tx counts 1 (unweighted). tau scaled by 0.9999 before arctanh. Single panel, poster-sized text, dpi=600.
gene tau_weighted = the same shared onebrain weight vector (IN_WEIGHTS_CSV) applied to the
gene's summed-isoform profile (as in compute_tx_info_content.py).
Single-panel poster version of plot_tau_deviation_histogram_weighted.py."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gls_tau import weighted_tau

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_WEIGHTS_CSV = sys.argv[3]
OUT_PNG = sys.argv[4]

SCALE = 0.9999
AXIS_LABEL_FONTSIZE = 20
TICK_LABEL_FONTSIZE = 16

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0
tx = pd.DataFrame({"transcript_id": gct["Name"][keep], "gene_id": gct["Description"][keep]})

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
w = pd.read_csv(IN_WEIGHTS_CSV).set_index("tissue").loc[tissues]["weight"].to_numpy()
gene_mat = gct.loc[keep, tissues].groupby(gct.loc[keep, "Description"]).sum()
info = pd.DataFrame({"gene_id": gene_mat.index,
                     "gene_tau_weighted": weighted_tau(gene_mat.to_numpy(), w)})

df = tx.merge(tx_tau, on="transcript_id").merge(info, on="gene_id")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted"])
n_tx = df.groupby("gene_id")["transcript_id"].transform("count")
df = df[n_tx >= 2].copy()
df["deviation"] = np.arctanh(SCALE * df["tx_tau_weighted"]) - np.arctanh(SCALE * df["gene_tau_weighted"])

fig, ax = plt.subplots(figsize=(10, 7))
ax.hist(df["deviation"], bins=100, color="steelblue", alpha=0.8)
ax.axvline(0, color="gray", linestyle="--", linewidth=2)
ax.set_xlabel("arctanh(tx tau_weighted) - arctanh(gene tau_weighted)", fontsize=AXIS_LABEL_FONTSIZE)
ax.set_ylabel("Transcript count", fontsize=AXIS_LABEL_FONTSIZE)
ax.tick_params(labelsize=TICK_LABEL_FONTSIZE)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)
print(f"n_tx={len(df)} n_genes={df['gene_id'].nunique()} median_dev={df['deviation'].median():.3f} -> {OUT_PNG}",
      file=sys.stderr)
