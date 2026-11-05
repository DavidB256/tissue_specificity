"""Scatter grid: rows = sources, cols = (i) tx total TPM (sum over tissues, log10), (ii) isoform fraction
(tx_total / gene_total). y = arctanh(0.9999*tx tau_w) - arctanh(0.9999*gene tau_w), multi-tx genes
(n_tx>=2), tx_total>0. gene tau_w = same onebrain weight vector applied to summed-isoform profile.
Usage: OUT_PNG label:gct:tx_tau_csv:weights_csv ..."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gls_tau import weighted_tau

OUT_PNG = sys.argv[1]
SOURCES = [a.split(":") for a in sys.argv[2:]]
SCALE = 0.9999


def load(gct_path, tx_tau_csv, weights_csv):
    gct = pd.read_csv(gct_path, sep="\t", skiprows=2)
    tissues = list(gct.columns[2:])
    tx_total = gct[tissues].sum(axis=1)
    keep = tx_total > 0
    tx = pd.DataFrame({"transcript_id": gct["Name"][keep], "gene_id": gct["Description"][keep],
                       "tx_total": tx_total[keep]})
    tx["isoform_fraction"] = tx["tx_total"] / tx.groupby("gene_id")["tx_total"].transform("sum")
    tx_tau = pd.read_csv(tx_tau_csv).rename(columns={"gene": "transcript_id", "tau_free_weighted": "tx_tau"})[
        ["transcript_id", "tx_tau"]]
    w = pd.read_csv(weights_csv).set_index("tissue").loc[tissues]["weight"].to_numpy()
    gene_mat = gct.loc[keep, tissues].groupby(gct.loc[keep, "Description"]).sum()
    gene = pd.DataFrame({"gene_id": gene_mat.index, "gene_tau": weighted_tau(gene_mat.to_numpy(), w)})
    df = tx.merge(tx_tau, on="transcript_id").merge(gene, on="gene_id").dropna(subset=["tx_tau", "gene_tau"])
    df = df[df.groupby("gene_id")["transcript_id"].transform("count") >= 2].copy()
    df["deviation"] = np.arctanh(SCALE * df["tx_tau"]) - np.arctanh(SCALE * df["gene_tau"])
    df["log10_tx_total"] = np.log10(df["tx_total"])
    return df


fig, axes = plt.subplots(len(SOURCES), 2, figsize=(14, 6 * len(SOURCES)), squeeze=False)
for r, (label, gct_path, tau_csv, w_csv) in enumerate(SOURCES):
    df = load(gct_path, tau_csv, w_csv)
    for c, (xcol, xlab) in enumerate([("log10_tx_total", "log10 transcript TPM (sum over tissues)"),
                                      ("isoform_fraction", "isoform fraction")]):
        ax = axes[r, c]
        ax.scatter(df[xcol], df["deviation"], s=6, alpha=0.1, color="#34688f", rasterized=True)
        ax.axhline(0, color="black", linestyle="--", linewidth=1.5)
        ax.set_xlabel(xlab)
        ax.set_ylabel("transcript tau_w - gene tau_w (arctanh)")
        ax.set_title(f"{label} (n={len(df)} tx)")
    print(f"{label}: n_tx={len(df)}", file=sys.stderr)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
