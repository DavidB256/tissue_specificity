"""Compare 8 tissue-clustering dendrograms built from every combination of
{pearson, spearman} x {TPM<1, TPM<10 floor} x {gtex_v11_tx_protein_coding,
chess_tx_protein_coding}. Each scenario: floor_tpm(mat, floor) -> tissue_corr_
expressed_only(method) -> average-linkage dendrogram (distance=1-corr) -> flat
clusters at distance threshold 0.1 (fcluster, same hierclust0.1 convention used
throughout this pipeline).

Dendrogram distance metric: 1 - cophenetic correlation. The cophenetic distance
between two leaves is just the merge height of their lowest common ancestor --
well-defined for any pair of leaves regardless of what other leaves are in the
tree, so for a pair of scenarios with DIFFERENT tissue sets (any gtex vs. chess
pair: 68 vs. 55 tissues, overlapping in 53 -- chess has 2 tissues gtex doesn't:
Cells_Transformed_fibroblasts, Cells_Leukemia_cell_line_CML) we restrict both
trees' cophenetic distances to the tissues present in both, then correlate
those two vectors. For a same-source pair the "common tissues"
set is just the full shared panel, reducing to the ordinary cophenetic
correlation. Also reports ARI on the flat hierclust0.1 cluster assignments,
restricted the same way, as a secondary comparison."""
import sys
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from sklearn.metrics import adjusted_rand_score
from gls_tau import tissue_corr_expressed_only, floor_tpm

IN_GTEX_GCT = sys.argv[1]
IN_CHESS_GCT = sys.argv[2]
OUT_CLUSTER_COUNTS_CSV = sys.argv[3]
OUT_DISTANCE_MATRIX_CSV = sys.argv[4]
OUT_DISTANCE_LONG_CSV = sys.argv[5]

DIST_THRESHOLD = 0.1
FLOORS = [1.0, 10.0]
METHODS = ["pearson", "spearman"]
SOURCES = {"gtex_v11_tx_protein_coding": IN_GTEX_GCT, "chess_tx_protein_coding": IN_CHESS_GCT}

gcts = {}
for source, path in SOURCES.items():
    df = pd.read_csv(path, sep="\t", skiprows=2)
    tissues = list(df.columns[2:])
    gcts[source] = (tissues, df[tissues].to_numpy(dtype=float))

common_tissues_all = sorted(set(gcts["gtex_v11_tx_protein_coding"][0]) & set(gcts["chess_tx_protein_coding"][0]))
print(f"tissues: gtex={len(gcts['gtex_v11_tx_protein_coding'][0])} "
      f"chess={len(gcts['chess_tx_protein_coding'][0])} common={len(common_tissues_all)}", file=sys.stderr)

scenarios = {}
rows = []
for source in SOURCES:
    tissues, raw_mat = gcts[source]
    for floor in FLOORS:
        mat = floor_tpm(raw_mat, floor=floor)
        for method in METHODS:
            name = f"{source}__floor{int(floor)}__{method}"
            R = tissue_corr_expressed_only(mat, method=method)
            dist = 1 - R
            np.fill_diagonal(dist, 0.0)
            dist = (dist + dist.T) / 2
            Z = linkage(squareform(dist, checks=False), method="average")
            cluster_id = fcluster(Z, t=DIST_THRESHOLD, criterion="distance")
            coph = squareform(cophenet(Z))
            n_clusters = len(set(cluster_id))
            scenarios[name] = {
                "tissues": tissues, "cophenetic": coph,
                "cluster_of": dict(zip(tissues, cluster_id)),
            }
            rows.append({
                "scenario": name, "source": source, "floor_threshold": floor,
                "corr_method": method, "n_tissues": len(tissues), "n_clusters": n_clusters,
            })
            print(f"{name}: n_tissues={len(tissues)} n_clusters={n_clusters}", file=sys.stderr)

pd.DataFrame(rows).to_csv(OUT_CLUSTER_COUNTS_CSV, index=False)

names = list(scenarios.keys())
dist_matrix = pd.DataFrame(np.zeros((len(names), len(names))), index=names, columns=names)
long_rows = []
for a, b in combinations(names, 2):
    sa, sb = scenarios[a], scenarios[b]
    common = sorted(set(sa["tissues"]) & set(sb["tissues"]))
    idx_a = {t: i for i, t in enumerate(sa["tissues"])}
    idx_b = {t: i for i, t in enumerate(sb["tissues"])}

    pairs = list(combinations(common, 2))
    coph_a = np.array([sa["cophenetic"][idx_a[t1], idx_a[t2]] for t1, t2 in pairs])
    coph_b = np.array([sb["cophenetic"][idx_b[t1], idx_b[t2]] for t1, t2 in pairs])
    r, _ = pearsonr(coph_a, coph_b) if len(pairs) >= 2 else (np.nan, np.nan)
    dist = 1 - r

    clusters_a = [sa["cluster_of"][t] for t in common]
    clusters_b = [sb["cluster_of"][t] for t in common]
    ari = adjusted_rand_score(clusters_a, clusters_b)

    dist_matrix.loc[a, b] = dist_matrix.loc[b, a] = dist
    long_rows.append({
        "scenario_a": a, "scenario_b": b, "n_common_tissues": len(common),
        "cophenetic_pearson_r": r, "distance_1_minus_r": dist, "ari": ari,
    })

dist_matrix.to_csv(OUT_DISTANCE_MATRIX_CSV)
pd.DataFrame(long_rows).sort_values("distance_1_minus_r").to_csv(OUT_DISTANCE_LONG_CSV, index=False)

print(f"-> {OUT_CLUSTER_COUNTS_CSV}, {OUT_DISTANCE_MATRIX_CSV}, {OUT_DISTANCE_LONG_CSV}", file=sys.stderr)
