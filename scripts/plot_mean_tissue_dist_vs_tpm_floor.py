"""Mean tissue-tissue distance (1 - Pearson r, off-diagonal mean) vs. TPM floor
(TPM<f -> 0, f=0..30; f=0 is a no-op, and is skipped for tx sources since their GCTs are already floored at TPM<1), one panel per source. Correlation as in
compute_tissue_correlation.py (log10(TPM+1), rows expressed in >=1 tissue after flooring).
95% bootstrap band: resample gene clusters (GCT Description column -- gene_id for tx
GCTs, so all isoforms of a gene move together; gene name for gene GCTs) with replacement,
same resamples at every floor (paired). Resampling is done via per-row multiplicity
weights (weighted Pearson == corrcoef on the duplicated matrix).
Usage: plot_mean_tissue_dist_vs_tpm_floor.py OUT_PNG OUT_CSV N_BOOT N_THREADS label:gct [label:gct ...]"""
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tau_utils import floor_tpm

OUT_PNG, OUT_CSV = sys.argv[1], sys.argv[2]
N_BOOT, N_THREADS = int(sys.argv[3]), int(sys.argv[4])
SOURCES = [a.split(":", 1) for a in sys.argv[5:]]
FLOORS = range(0, 31)
SEED = 0


def mean_offdiag_dist(X, w):
    """X: rows x tissues log10(TPM+1); w: row multiplicities. Weighted-Pearson mean of 1 - r."""
    sw = w.sum()
    mu = (w @ X) / sw
    C = (X * w[:, None]).T @ X / sw - np.outer(mu, mu)
    sd = np.sqrt(np.diag(C))
    R = C / np.outer(sd, sd)
    off = ~np.eye(len(R), dtype=bool)
    return float((1 - R[off]).mean())


rng = np.random.default_rng(SEED)
rows = []
for label, gct in SOURCES:
    df = pd.read_csv(gct, sep="\t", skiprows=2)
    raw = df[list(df.columns[2:])].values.astype(float)
    cluster = pd.factorize(df["Description"])[0]
    n_clu = cluster.max() + 1
    counts = np.stack([np.bincount(rng.integers(0, n_clu, n_clu), minlength=n_clu) for _ in range(N_BOOT)])
    for f in FLOORS:
        if f == 0 and "tx" in label:
            continue
        floored = floor_tpm(raw, floor=float(f))
        expressed = floored.max(axis=1) > 0
        X = np.log10(floored[expressed] + 1.0)
        cl = cluster[expressed]
        point = mean_offdiag_dist(X, np.ones(len(X)))
        with ThreadPoolExecutor(N_THREADS) as ex:
            boot = np.array(list(ex.map(lambda c: mean_offdiag_dist(X, c[cl].astype(float)), counts)))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rows.append((label, f, point, lo, hi, boot.std()))
        print(rows[-1], file=sys.stderr)

res = pd.DataFrame(rows, columns=["source", "tpm_floor", "mean_dist", "ci_lo", "ci_hi", "boot_sd"])
res.to_csv(OUT_CSV, index=False)

fig, axes = plt.subplots(1, len(SOURCES), figsize=(5 * len(SOURCES), 4.5), sharey=True)
for ax, (label, _) in zip(np.atleast_1d(axes), SOURCES):
    sub = res[res["source"] == label]
    ax.fill_between(sub["tpm_floor"], sub["ci_lo"], sub["ci_hi"], alpha=0.3)
    ax.plot(sub["tpm_floor"], sub["mean_dist"], marker="o")
    ax.set_title(label)
    ax.set_xlabel("TPM floor (TPM < x -> 0)")
    ax.set_xticks(list(range(0, 31, 5)))
np.atleast_1d(axes)[0].set_ylabel("mean tissue distance (1 - r)")
fig.suptitle(f"95% bootstrap band, {N_BOOT} gene-cluster resamples", fontsize=10)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
