"""Plot Exponential+2Gaussian fits (best/median/worst by KS) and BIC comparison vs GMM k=3."""
import sys
import pandas as pd
import numpy as np
from scipy.stats import expon, norm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_PARQUET = sys.argv[1]
IN_FITS = sys.argv[2]
IN_GMM_FITS = sys.argv[3]
OUT_GRID_PNG = sys.argv[4]
OUT_COMPARE_PNG = sys.argv[5]

long_df = pd.read_parquet(IN_PARQUET)
fits = pd.read_csv(IN_FITS)
gmm_fits = pd.read_csv(IN_GMM_FITS)

fits_sorted = fits.sort_values("ks_stat").reset_index(drop=True)
n = len(fits_sorted)
picks = {
    "best fit": fits_sorted.iloc[0],
    "median fit": fits_sorted.iloc[n // 2],
    "worst fit": fits_sorted.iloc[-1],
}

fig, axes = plt.subplots(1, 3, figsize=(15, 4))
for ax, (label, row) in zip(axes, picks.items()):
    x = long_df.loc[long_df.tissue == row.tissue, "log10_tpm"].values
    ax.hist(x, bins=40, density=True, alpha=0.5, color="steelblue", label="data")
    xs = np.linspace(x.min(), x.max(), 500)
    comp_expon = row.w_expon * expon.pdf(xs, scale=1 / row.rate)
    comp_n1 = row.w_norm1 * norm.pdf(xs, loc=row.mean1, scale=row.std1)
    comp_n2 = row.w_norm2 * norm.pdf(xs, loc=row.mean2, scale=row.std2)
    total = comp_expon + comp_n1 + comp_n2
    ax.plot(xs, comp_expon, "--", lw=1, color="green", label="expon")
    ax.plot(xs, comp_n1, "--", lw=1, color="gray", label="norm1")
    ax.plot(xs, comp_n2, "--", lw=1, color="gray", label="norm2")
    ax.plot(xs, total, color="darkred", lw=2, label="mixture")
    ax.set_title(f"{label}: {row.tissue}\nks={row.ks_stat:.4f} p={row.ks_pval:.2g} bic={row.bic:.0f}", fontsize=9)
    ax.set_xlabel("log10(TPM)")
    ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT_GRID_PNG, dpi=300)
plt.close(fig)

# BIC comparison: expgauss vs best GMM (k=1,2,3) per tissue
gmm_best = gmm_fits.loc[gmm_fits.groupby("tissue")["bic"].idxmin(), ["tissue", "bic", "ks_stat"]]
gmm_best = gmm_best.rename(columns={"bic": "bic_gmm", "ks_stat": "ks_gmm"})
merged = fits.merge(gmm_best, on="tissue").rename(columns={"bic": "bic_expgauss", "ks_stat": "ks_expgauss"})

fig2, axes2 = plt.subplots(1, 2, figsize=(10, 4.5))
lims_bic = [merged[["bic_gmm", "bic_expgauss"]].min().min(), merged[["bic_gmm", "bic_expgauss"]].max().max()]
axes2[0].scatter(merged.bic_gmm, merged.bic_expgauss, s=15, alpha=0.7)
axes2[0].plot(lims_bic, lims_bic, "k--", lw=1)
axes2[0].set_xlabel("BIC: best GMM (k=1,2,3)")
axes2[0].set_ylabel("BIC: Exponential + 2 Gaussian")
axes2[0].set_title("BIC per tissue (below diagonal = expgauss better)")

lims_ks = [0, merged[["ks_gmm", "ks_expgauss"]].max().max()]
axes2[1].scatter(merged.ks_gmm, merged.ks_expgauss, s=15, alpha=0.7)
axes2[1].plot(lims_ks, lims_ks, "k--", lw=1)
axes2[1].set_xlabel("KS stat: best GMM (k=1,2,3)")
axes2[1].set_ylabel("KS stat: Exponential + 2 Gaussian")
axes2[1].set_title("KS stat per tissue (below diagonal = expgauss better)")
fig2.tight_layout()
fig2.savefig(OUT_COMPARE_PNG, dpi=300)
plt.close(fig2)

n_better_bic = (merged.bic_expgauss < merged.bic_gmm).sum()
n_better_ks = (merged.ks_expgauss < merged.ks_gmm).sum()
n_nonsig = (merged.ks_pval > 0.05).sum()
print(f"expgauss better BIC: {n_better_bic}/{len(merged)}", file=sys.stderr)
print(f"expgauss better KS: {n_better_ks}/{len(merged)}", file=sys.stderr)
print(f"expgauss KS non-significant (p>0.05): {n_nonsig}/{len(merged)}", file=sys.stderr)
