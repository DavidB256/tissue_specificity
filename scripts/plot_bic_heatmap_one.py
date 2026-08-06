"""Per-tissue fit-quality landscape: two 20x20 heatmaps over disjoint (rank(log10TPM),
rank(tau)) bins -- a genuine 2D histogram, each cell an equal-width (not equal-count,
since u,v are correlated) rank interval on each axis. For each cell, using only the
points that fall in that cell:
  exp LL/n  = mean log-likelihood per point of fitting a lone exponential to the cell
  beta LL/n = mean log-likelihood per point of fitting a lone (bounded) beta to the cell
LL/n (not BIC) so cells are comparable regardless of how many points they contain --
this shows which regions of the rank plane locally look exponential- vs beta-shaped.
Cells with <10 points are masked."""
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

bin_u = np.clip((u * N_BINS).astype(int), 0, N_BINS - 1)
bin_v = np.clip((v * N_BINS).astype(int), 0, N_BINS - 1)

exp_ll = np.full((N_BINS, N_BINS), np.nan)
beta_ll = np.full((N_BINS, N_BINS), np.nan)

for i in range(N_BINS):
    for j in range(N_BINS):
        cell_mask = (bin_u == i) & (bin_v == j)
        n_cell = int(cell_mask.sum())
        if n_cell < MIN_N:
            continue
        x_cell = x[cell_mask]

        rate = 1.0 / x_cell.mean()
        exp_ll[j, i] = logpdf(x_cell, "expon", {"rate": rate}).mean()

        beta_fit = fit_mixture(x_cell, ["beta"], n_init=2)
        if beta_fit["loglik"] > -1e6:  # else: optimizer hit the non-finite-likelihood fallback penalty
            beta_ll[j, i] = beta_fit["loglik"] / n_cell

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, grid, title in [(axes[0], exp_ll, "exp: mean log-likelihood / point"),
                         (axes[1], beta_ll, "beta: mean log-likelihood / point")]:
    im = ax.imshow(grid, origin="lower", extent=[0, 1, 0, 1], aspect="auto", cmap="viridis")
    fig.colorbar(im, ax=ax, label="LL/n")
    ax.set_xlabel("rank(log10TPM) bin")
    ax.set_ylabel("rank(tau) bin")
    ax.set_title(title, fontsize=10)

fig.suptitle(f"{TISSUE}: local fit quality (LL/n) over disjoint (log10TPM rank, tau rank) bins, {N_BINS}x{N_BINS}")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
