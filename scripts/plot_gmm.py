"""Plot histogram + best-k GMM fit for best/median/worst-KS tissues, and BIC-vs-k summary."""
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_PARQUET = sys.argv[1]
IN_FITS = sys.argv[2]
OUT_GRID_PNG = sys.argv[3]
OUT_BIC_PNG = sys.argv[4]

long_df = pd.read_parquet(IN_PARQUET)
fits = pd.read_csv(IN_FITS)

# pick best k per tissue by BIC
best = fits.loc[fits.groupby("tissue")["bic"].idxmin()].reset_index(drop=True)

# rank tissues by ks_stat of their best-k fit: best / median / worst
best_sorted = best.sort_values("ks_stat").reset_index(drop=True)
n = len(best_sorted)
picks = {
    "best fit": best_sorted.iloc[0],
    "median fit": best_sorted.iloc[n // 2],
    "worst fit": best_sorted.iloc[-1],
}

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, (label, row) in zip(axes, picks.items()):
    x = long_df.loc[long_df.tissue == row.tissue, "log10_tpm"].values
    ax.hist(x, bins=40, density=True, alpha=0.5, color="steelblue", label="data")
    means = [float(v) for v in row["means"].split(";")]
    stds = [float(v) for v in row["stds"].split(";")]
    weights = [float(v) for v in row["weights"].split(";")]
    xs = np.linspace(x.min(), x.max(), 500)
    total = np.zeros_like(xs)
    for m, s, w in zip(means, stds, weights):
        comp = w * (1 / (s * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((xs - m) / s) ** 2)
        total += comp
        ax.plot(xs, comp, "--", lw=1, color="gray")
    ax.plot(xs, total, color="darkred", lw=2, label=f"GMM k={int(row.k)}")
    ax.set_title(f"{label}: {row.tissue}\nks={row.ks_stat:.3f} bic={row.bic:.0f}", fontsize=9)
    ax.set_xlabel("log10(TPM)")
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_GRID_PNG, dpi=300)
plt.close(fig)

# BIC vs k, averaged across tissues, to see if mixture components are ever justified
fig2, ax2 = plt.subplots(figsize=(5, 4))
mean_bic = fits.groupby("k")["bic"].mean()
ax2.plot(mean_bic.index, mean_bic.values, "o-")
ax2.set_xlabel("k (number of components)")
ax2.set_ylabel("mean BIC across tissues")
ax2.set_title("GMM: mean BIC vs k")
ax2.set_xticks(mean_bic.index)
fig2.tight_layout()
fig2.savefig(OUT_BIC_PNG, dpi=300)
plt.close(fig2)

print(best["k"].value_counts(), file=sys.stderr)
