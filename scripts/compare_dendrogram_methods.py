"""How similar are two hierarchical tissue clusterings built from different
correlation metrics (e.g. pearson vs spearman)? Both linkages were built from
the same GCT / tissue column order, so leaf index i is the same tissue in
both -- cophenetic distances are directly comparable without re-matching
leaves. Clusters CSVs are matched by tissue name (join), since cluster_id
numbering is arbitrary per-tree."""
import sys
import pandas as pd
from scipy.cluster.hierarchy import cophenet
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import adjusted_rand_score, rand_score

IN_LINKAGE_A = sys.argv[1]
IN_LINKAGE_B = sys.argv[2]
IN_CLUSTERS_A = sys.argv[3]
IN_CLUSTERS_B = sys.argv[4]
OUT_CSV = sys.argv[5]

Z_a = pd.read_csv(IN_LINKAGE_A)[["cluster1", "cluster2", "distance", "n_obs"]].values
Z_b = pd.read_csv(IN_LINKAGE_B)[["cluster1", "cluster2", "distance", "n_obs"]].values

coph_a = cophenet(Z_a)
coph_b = cophenet(Z_b)
coph_pearson_r, coph_pearson_p = pearsonr(coph_a, coph_b)
coph_spearman_rho, coph_spearman_p = spearmanr(coph_a, coph_b)

clusters_a = pd.read_csv(IN_CLUSTERS_A)[["tissue", "cluster_id"]].rename(columns={"cluster_id": "cluster_a"})
clusters_b = pd.read_csv(IN_CLUSTERS_B)[["tissue", "cluster_id"]].rename(columns={"cluster_id": "cluster_b"})
merged = clusters_a.merge(clusters_b, on="tissue", how="inner")

ari = adjusted_rand_score(merged["cluster_a"], merged["cluster_b"])
ri = rand_score(merged["cluster_a"], merged["cluster_b"])

out = pd.DataFrame([{
    "n_tissues": len(merged),
    "n_clusters_a": merged["cluster_a"].nunique(),
    "n_clusters_b": merged["cluster_b"].nunique(),
    "cophenetic_pearson_r": coph_pearson_r,
    "cophenetic_pearson_p": coph_pearson_p,
    "cophenetic_spearman_rho": coph_spearman_rho,
    "cophenetic_spearman_p": coph_spearman_p,
    "adjusted_rand_index": ari,
    "rand_index": ri,
}])
out.to_csv(OUT_CSV, index=False)

print(f"n_tissues={len(merged)} n_clusters_a={merged['cluster_a'].nunique()} "
      f"n_clusters_b={merged['cluster_b'].nunique()} cophenetic_pearson_r={coph_pearson_r:.4f} "
      f"ARI={ari:.4f} RI={ri:.4f} -> {OUT_CSV}", file=sys.stderr)
