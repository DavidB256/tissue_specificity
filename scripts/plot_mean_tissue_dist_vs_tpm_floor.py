"""Mean tissue-tissue distance (1 - Pearson r, off-diagonal mean) vs. TPM floor
(TPM<f -> 0, f=0..20; f=0 is a no-op), one panel per source. Correlation as in
compute_tissue_correlation.py (log10(TPM+1), genes expressed in >=1 tissue after flooring).
Usage: plot_mean_tissue_dist_vs_tpm_floor.py OUT_PNG OUT_CSV label:gct [label:gct ...]"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gls_tau import tissue_corr_expressed_only, floor_tpm

OUT_PNG, OUT_CSV = sys.argv[1], sys.argv[2]
SOURCES = [a.split(":", 1) for a in sys.argv[3:]]
FLOORS = range(0, 21)

rows = []
for label, gct in SOURCES:
    df = pd.read_csv(gct, sep="\t", skiprows=2)
    raw = df[list(df.columns[2:])].values.astype(float)
    for f in FLOORS:
        R = tissue_corr_expressed_only(floor_tpm(raw, floor=float(f)))
        off = ~np.eye(len(R), dtype=bool)
        rows.append((label, f, float((1 - R[off]).mean())))
        print(rows[-1], file=sys.stderr)

res = pd.DataFrame(rows, columns=["source", "tpm_floor", "mean_dist"])
res.to_csv(OUT_CSV, index=False)

fig, axes = plt.subplots(1, len(SOURCES), figsize=(5 * len(SOURCES), 4.5), sharey=True)
for ax, (label, _) in zip(np.atleast_1d(axes), SOURCES):
    sub = res[res["source"] == label]
    ax.plot(sub["tpm_floor"], sub["mean_dist"], marker="o")
    ax.set_title(label)
    ax.set_xlabel("TPM floor (TPM < x -> 0)")
    ax.set_xticks(list(range(0, 21, 2)))
np.atleast_1d(axes)[0].set_ylabel("mean tissue distance (1 - r)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
