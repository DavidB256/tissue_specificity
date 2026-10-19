"""Mantel test: does each source's transcriptomic tissue-tissue cophenetic
distance (from its average-linkage dendrogram, tissue_correlation_linkage.csv)
agree with UBERON anatomical distance (1-batet, uberon_tissue_distance.csv)?
Recommended over forcing UBERON into a binary tree + Robinson-Foulds (see
uberon_dag_tree_comparison_ideas.md option 5): skips building any tree on the
UBERON side (avoiding the ultrametric distortion average-linkage would
introduce before an RF comparison), and tests the two distance structures'
correlation directly.

For each source: restrict both matrices to their common tissues, correlate
the upper-triangle distance vectors (Pearson r = "Mantel statistic"), then
permutation-test significance by randomly relabeling (permuting) the tissue
order of one matrix N times and recomputing r each time -- p = fraction of
permuted |r| >= observed |r| (two-sided), +1/+1 smoothing. This correctly
accounts for the non-independence of distance-matrix entries, which a naive
pairwise-correlation p-value would ignore."""
import sys
import pickle
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr

IN_UBERON_DISTANCE_CSV = sys.argv[1]
OUT_RESULTS_CSV = sys.argv[2]
SOURCES = sys.argv[3:]  # "name:tissue_correlation_csv:linkage_csv" triples

N_PERM = 9999
RNG = np.random.default_rng(0)

uberon = pd.read_csv(IN_UBERON_DISTANCE_CSV, index_col=0)

rows = []
mantel_data = {}
for spec in SOURCES:
    name, corr_csv, linkage_csv = spec.split(":")
    corr = pd.read_csv(corr_csv, index_col=0)
    tissues = list(corr.index)
    Z = pd.read_csv(linkage_csv)[["cluster1", "cluster2", "distance", "n_obs"]].to_numpy()
    coph = pd.DataFrame(squareform(cophenet(Z)), index=tissues, columns=tissues)

    common = sorted(set(tissues) & set(uberon.index))
    coph_sub = coph.loc[common, common].to_numpy()
    uberon_sub = uberon.loc[common, common].to_numpy()
    iu = np.triu_indices(len(common), k=1)
    coph_vec = coph_sub[iu]
    uberon_vec = uberon_sub[iu]

    r_obs, _ = pearsonr(uberon_vec, coph_vec)

    n = len(common)
    count_extreme = 0
    for _ in range(N_PERM):
        perm = RNG.permutation(n)
        coph_perm = coph_sub[np.ix_(perm, perm)][iu]
        r_perm, _ = pearsonr(uberon_vec, coph_perm)
        if abs(r_perm) >= abs(r_obs):
            count_extreme += 1
    p_value = (count_extreme + 1) / (N_PERM + 1)

    mantel_data[name] = {"uberon_dist": uberon_vec, "coph_dist": coph_vec}
    rows.append({"source": name, "n_common_tissues": n, "mantel_r": r_obs, "p_value": p_value})
    print(f"{name}: n_common={n} mantel_r={r_obs:.4f} p={p_value:.4f}", file=sys.stderr)

pd.DataFrame(rows).to_csv(OUT_RESULTS_CSV, index=False)

with open(OUT_RESULTS_CSV.replace(".csv", "_vectors.pkl"), "wb") as f:
    pickle.dump(mantel_data, f)

print(f"-> {OUT_RESULTS_CSV}", file=sys.stderr)
