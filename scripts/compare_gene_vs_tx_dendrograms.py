"""Compare dendrograms across the gene-centric source (gtex_protein_coding)
and the transcript-centric sources (gtex_v11_tx_protein_coding, chess_tx_
protein_coding), at both TPM<1 and TPM<10 flooring, so the floor's effect can
be compared alongside the gene-vs-tx effect. floor_tpm(mat, floor) -> tissue_
corr_expressed_only(method="pearson") (which, after flooring, is automatically
restricted to genes/tx reaching >=floor TPM in at least one tissue -- the
"max(axis=1) > 0" filter already baked into tissue_corr_expressed_only
operates on the POST-floor matrix, so nothing extra is needed to satisfy that
exclusion). Pearson only, per instruction -- no spearman. Flat clusters at
distance threshold 0.1, same hierclust0.1 convention used throughout this
pipeline.

baron_protein_coding is deliberately excluded here (it shares no leaf names
with the other sources -- GTEx tissue names vs. Baron pancreas cell-type
names -- so every comparison involving it was NaN; see
plot_baron_umap_thresholds.py / the earlier gene_vs_tx_dendrograms.png for
Baron-specific floor-threshold analysis instead)."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
from gls_tau import tissue_corr_expressed_only, floor_tpm
from tree_compare_utils import make_scenario, compare_all_scenarios

IN_GTEX_PC_GCT = sys.argv[1]
IN_GTEX_TX_PC_GCT = sys.argv[2]
IN_CHESS_TX_PC_GCT = sys.argv[3]
OUT_CLUSTER_COUNTS_CSV = sys.argv[4]
OUT_DISTANCE_MATRIX_CSV = sys.argv[5]
OUT_RF_MATRIX_CSV = sys.argv[6]
OUT_DISTANCE_LONG_CSV = sys.argv[7]

FLOORS = [1.0, 10.0]
METHOD = "pearson"
SOURCES = {
    "gtex_protein_coding": ("gene_centric", IN_GTEX_PC_GCT),
    "gtex_v11_tx_protein_coding": ("tx_centric", IN_GTEX_TX_PC_GCT),
    "chess_tx_protein_coding": ("tx_centric", IN_CHESS_TX_PC_GCT),
}

scenarios = {}
rows = []
for source, (kind, path) in SOURCES.items():
    df = pd.read_csv(path, sep="\t", skiprows=2)
    tissues = list(df.columns[2:])
    raw_mat = df[tissues].to_numpy(dtype=float)
    for floor in FLOORS:
        mat = floor_tpm(raw_mat, floor=floor)
        n_pass = int((mat.max(axis=1) > 0).sum())

        R = tissue_corr_expressed_only(mat, method=METHOD)
        dist = 1 - R
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2
        Z = linkage(squareform(dist, checks=False), method="average")
        name = f"{source}__floor{int(floor)}"
        scenarios[name] = make_scenario(Z, tissues)
        rows.append({
            "scenario": name, "source": source, "kind": kind, "floor_threshold": floor,
            "n_tissues": len(tissues), "n_genes_or_tx_ge_floor_tpm": n_pass,
            "n_clusters": scenarios[name]["n_clusters"],
        })
        print(f"{name} ({kind}): n_tissues={len(tissues)} n_genes_or_tx>={int(floor)}TPM_in_1+tissue={n_pass} "
              f"n_clusters={scenarios[name]['n_clusters']}", file=sys.stderr)

pd.DataFrame(rows).to_csv(OUT_CLUSTER_COUNTS_CSV, index=False)

dist_matrix, rf_matrix, long_df = compare_all_scenarios(scenarios)
dist_matrix.to_csv(OUT_DISTANCE_MATRIX_CSV)
rf_matrix.to_csv(OUT_RF_MATRIX_CSV)
long_df.sort_values("distance_1_minus_r").to_csv(OUT_DISTANCE_LONG_CSV, index=False)

print(f"-> {OUT_CLUSTER_COUNTS_CSV}, {OUT_DISTANCE_MATRIX_CSV}, {OUT_RF_MATRIX_CSV}, {OUT_DISTANCE_LONG_CSV}",
      file=sys.stderr)
