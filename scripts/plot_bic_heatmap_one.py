"""Per-tissue BIC-landscape diagnostic: two 20x20 heatmaps over (rank(log10TPM),
rank(tau)) threshold pairs (u_cut, v_cut). For each grid cell:
  exp BIC  = BIC of fitting a lone exponential to the corner box
             {rank_tpm <= u_cut AND rank_tau <= v_cut}
  beta BIC = BIC of fitting a lone (bounded) beta to the complementary corner box
             {rank_tpm >= u_cut AND rank_tau >= v_cut}
Each BIC is for that single distribution alone (k=1 for exp, k=4 for beta), not a
joint two-component score -- these are meant to show whether BIC improves smoothly
as the corner sweeps the rank plane, independent of what the earlier line-search
(angle, offset) grid search actually picked. Cells with <10 points are masked."""
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mixture_core import fit_mixture, logpdf

warnings.filterwarnings("ignore", category=RuntimeWarning)

TISSUE = sys.argv[1]
IN_PARQUET = sys.argv[2]
IN_PARTITION_LABELS = sys.argv[3]
OUT_PNG = sys.argv[4]

N_BINS = 20
MIN_N = 10

long_df = pd.read_parquet(IN_PARQUET)
x_all = long_df.loc[long_df.tissue == TISSUE, ["gene", "log10_tpm"]]
labels = pd.read_csv(IN_PARTITION_LABELS)[["gene", "rank_tpm", "rank_tau"]]
merged = x_all.merge(labels, on="gene", how="inner")

x = merged["log10_tpm"].values
u = merged["rank_tpm"].values
v = merged["rank_tau"].values

thresholds = np.linspace(0.05, 1.0, N_BINS)
exp_bic = np.full((N_BINS, N_BINS), np.nan)
beta_bic = np.full((N_BINS, N_BINS), np.nan)

for i, u_cut in enumerate(thresholds):
    for j, v_cut in enumerate(thresholds):
        exp_mask = (u <= u_cut) & (v <= v_cut)
        beta_mask = (u >= u_cut) & (v >= v_cut)

        n_exp = int(exp_mask.sum())
        if n_exp >= MIN_N:
            x_exp = x[exp_mask]
            rate = 1.0 / x_exp.mean()
            ll_exp = logpdf(x_exp, "expon", {"rate": rate}).sum()
            exp_bic[j, i] = 1 * np.log(n_exp) - 2 * ll_exp

        n_beta = int(beta_mask.sum())
        if n_beta >= MIN_N:
            x_beta = x[beta_mask]
            beta_fit = fit_mixture(x_beta, ["beta"], n_init=2)
            if beta_fit["loglik"] > -1e6:  # else: optimizer hit the non-finite-likelihood fallback penalty
                beta_bic[j, i] = 4 * np.log(n_beta) - 2 * beta_fit["loglik"]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, grid, title in [(axes[0], exp_bic, "exp BIC (lower-left corner box)"),
                         (axes[1], beta_bic, "beta BIC (upper-right corner box)")]:
    im = ax.imshow(grid, origin="lower", extent=[0.05, 1.0, 0.05, 1.0], aspect="auto", cmap="viridis_r")
    fig.colorbar(im, ax=ax, label="BIC")
    ax.set_xlabel("u_cut = rank(log10TPM) threshold")
    ax.set_ylabel("v_cut = rank(tau) threshold")
    ax.set_title(title, fontsize=10)

fig.suptitle(f"{TISSUE}: BIC landscape over (log10TPM rank, tau rank) corner thresholds, {N_BINS}x{N_BINS}")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
