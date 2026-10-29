"""Gene-cluster bootstrap 95% CIs for compute_gene_tx_uberon_mantel.py's Mantel r
(GTEx gene- and tx-level sources x TPM<1/TPM<10 floors). Each resample redraws gene
clusters (GCT Description = gene_id, so a gene's transcripts move together) with
replacement; row multiplicities enter a weighted Pearson (== corrcoef on the duplicated
matrix) of log10(TPM+1) across tissues, then 1-r -> average-linkage -> cophenetic ->
Mantel r vs UBERON on the common tissues. Same resamples reused across floors.
Usage: bootstrap_gene_tx_uberon_mantel.py UBERON_CSV OUT_CSV N_BOOT N_THREADS label:source:gct ..."""
import os
for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[v] = "1"  # parallelism is across resamples
import sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform
from gls_tau import floor_tpm

IN_UBERON_DISTANCE_CSV = sys.argv[1]
OUT_CSV = sys.argv[2]
N_BOOT, N_THREADS = int(sys.argv[3]), int(sys.argv[4])
SOURCES = [a.split(":", 2) for a in sys.argv[5:]]
FLOORS = [1.0, 10.0]
SEED = 0


def weighted_corr(X, w):
    sw = w.sum()
    mu = (w @ X) / sw
    C = (X * w[:, None]).T @ X / sw - np.outer(mu, mu)
    sd = np.sqrt(np.diag(C))
    return C / np.outer(sd, sd)


def mantel_r(R, idx, uberon_vec, iu):
    dist = 1 - R
    np.fill_diagonal(dist, 0.0)
    dist = (dist + dist.T) / 2
    coph = squareform(cophenet(linkage(squareform(dist, checks=False), method="average")))
    return np.corrcoef(uberon_vec, coph[np.ix_(idx, idx)][iu])[0, 1]


uberon = pd.read_csv(IN_UBERON_DISTANCE_CSV, index_col=0)
rng = np.random.default_rng(SEED)
rows = []
for label, source, gct in SOURCES:
    df = pd.read_csv(gct, sep="\t", skiprows=2)
    tissues = list(df.columns[2:])
    raw = df[tissues].to_numpy(dtype=float)
    cluster = pd.factorize(df["Description"])[0]
    n_clu = cluster.max() + 1
    counts = np.stack([np.bincount(rng.integers(0, n_clu, n_clu), minlength=n_clu) for _ in range(N_BOOT)])
    common = sorted(set(tissues) & set(uberon.index))
    idx = np.array([tissues.index(t) for t in common])
    iu = np.triu_indices(len(common), k=1)
    uberon_vec = uberon.loc[common, common].to_numpy()[iu]
    for floor in FLOORS:
        floored = floor_tpm(raw, floor=floor)
        expressed = floored.max(axis=1) > 0
        X = np.log10(floored[expressed] + 1.0)
        cl = cluster[expressed]
        point = mantel_r(weighted_corr(X, np.ones(len(X))), idx, uberon_vec, iu)
        with ThreadPoolExecutor(N_THREADS) as ex:
            boot = np.array(list(ex.map(
                lambda c: mantel_r(weighted_corr(X, c[cl].astype(float)), idx, uberon_vec, iu), counts)))
        lo, hi = np.percentile(boot, [2.5, 97.5])
        rows.append({"scenario": f"{source}__floor{int(floor)}", "source": source, "floor_threshold": floor,
                     "mantel_r": point, "ci_lo": lo, "ci_hi": hi, "boot_sd": boot.std(), "n_boot": N_BOOT})
        print(rows[-1], file=sys.stderr)

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f"-> {OUT_CSV}", file=sys.stderr)
