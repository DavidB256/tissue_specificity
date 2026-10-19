"""Compare 13 tissue-clustering dendrograms: 12 built from every combination of
{pearson, spearman} x {none, TPM<1, TPM<10 floor} x {gtex_v11_tx_protein_coding,
chess_tx_protein_coding}, plus 1 from UBERON ontology distance (gtex tissues
only -- see compute_uberon_tissue_distance.py). "none" (floor=0) never rounds
any value down to 0 -- floor_tpm(mat, floor=0) is a no-op since TPM is always
>=0. Each correlation scenario: floor_tpm(mat, floor) -> tissue_corr_
expressed_only(method) -> average-linkage dendrogram (distance=1-corr) -> flat
clusters at distance threshold 0.1 (fcluster, same hierclust0.1 convention
used throughout this pipeline) -- note this threshold was calibrated for the
1-corr distance scale and is NOT meaningful for the UBERON scenario's batet-
based distances (expect its n_clusters to land near 68, i.e. mostly
singletons); the cophenetic/RF comparisons below don't depend on this cut.

Two distinct notions of "different dendrograms" are measured, restricted to
common tissues the same way for cross-source pairs (any gtex vs. chess pair:
68 vs. 55 tissues, overlapping in 53 -- chess has 2 tissues gtex doesn't:
Cells_Transformed_fibroblasts, Cells_Leukemia_cell_line_CML):

1. Cophenetic distance (distance_1_minus_r column + the distance-matrix CSV):
   1 - Pearson correlation between the two trees' cophenetic distances (merge
   height of each pair of leaves' lowest common ancestor) over common-tissue
   pairs. Height-sensitive: trees with the same branching order but different
   merge heights still show up as "different" here.
2. Robinson-Foulds topology distance (rf_jaccard_distance column + the
   rf-matrix CSV): every internal node of a tree induces a bipartition of the
   leaves (its descendants vs. everyone else); collect all NON-trivial ones
   (2..n-2 leaves on the smaller side) per tree, restrict each to the common
   tissues (re-canonicalized, since restricting can shrink/trivialize a split),
   and compare the two trees' bipartition SETS via Jaccard distance
   (1 - |intersection|/|union|). This ignores merge heights entirely and asks
   only "do the same groupings of tissues appear in both trees" -- pure
   topology, complementary to the height-sensitive cophenetic metric.

Also reports ARI on the flat hierclust0.1 cluster assignments, restricted the
same way, as a third, single-cut-only comparison."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from gls_tau import tissue_corr_expressed_only, floor_tpm
from tree_compare_utils import make_scenario, compare_all_scenarios

IN_GTEX_GCT = sys.argv[1]
IN_CHESS_GCT = sys.argv[2]
IN_UBERON_LINKAGE_CSV = sys.argv[3]
IN_UBERON_MAPPING_CSV = sys.argv[4]
OUT_CLUSTER_COUNTS_CSV = sys.argv[5]
OUT_DISTANCE_MATRIX_CSV = sys.argv[6]
OUT_RF_MATRIX_CSV = sys.argv[7]
OUT_DISTANCE_LONG_CSV = sys.argv[8]

DIST_THRESHOLD = 0.1
FLOORS = [0.0, 1.0, 10.0]
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
            scenarios[name] = make_scenario(Z, tissues)
            rows.append({
                "scenario": name, "source": source, "floor_threshold": floor,
                "corr_method": method, "n_tissues": len(tissues),
                "n_clusters": scenarios[name]["n_clusters"],
            })
            print(f"{name}: n_tissues={len(tissues)} n_clusters={scenarios[name]['n_clusters']}", file=sys.stderr)

uberon_tissues = list(pd.read_csv(IN_UBERON_MAPPING_CSV)["tissue"])
uberon_Z = pd.read_csv(IN_UBERON_LINKAGE_CSV)[["cluster1", "cluster2", "distance", "n_obs"]].values
scenarios["gtex_uberon_tree"] = make_scenario(uberon_Z, uberon_tissues)
rows.append({
    "scenario": "gtex_uberon_tree", "source": "gtex_uberon_tree", "floor_threshold": np.nan,
    "corr_method": "n/a", "n_tissues": len(uberon_tissues),
    "n_clusters": scenarios["gtex_uberon_tree"]["n_clusters"],
})
print(f"gtex_uberon_tree: n_tissues={len(uberon_tissues)} "
      f"n_clusters={scenarios['gtex_uberon_tree']['n_clusters']}", file=sys.stderr)

pd.DataFrame(rows).to_csv(OUT_CLUSTER_COUNTS_CSV, index=False)

dist_matrix, rf_matrix, long_df = compare_all_scenarios(scenarios)
dist_matrix.to_csv(OUT_DISTANCE_MATRIX_CSV)
rf_matrix.to_csv(OUT_RF_MATRIX_CSV)
long_df.sort_values("distance_1_minus_r").to_csv(OUT_DISTANCE_LONG_CSV, index=False)

print(f"-> {OUT_CLUSTER_COUNTS_CSV}, {OUT_DISTANCE_MATRIX_CSV}, {OUT_RF_MATRIX_CSV}, {OUT_DISTANCE_LONG_CSV}",
      file=sys.stderr)
