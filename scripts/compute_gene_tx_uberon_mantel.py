"""Mantel test (see compute_uberon_mantel.py for the method) of each of the 6
gene_vs_tx_bars.png scenarios (3 sources: gtex_protein_coding [gene-level],
gtex_v11_tx_protein_coding/chess_tx_protein_coding [tx-level] x 2 TPM floors:
<1, <10) against UBERON anatomical (1-batet) distance. For each scenario:
floor_tpm -> tissue_corr_expressed_only(pearson) -> average-linkage dendrogram
-> cophenetic distance matrix -> restrict to tissues common with UBERON ->
Mantel r (Pearson r between the two distance matrices' upper-triangle
entries) with a 9999-permutation two-sided p-value. Adds a "stars" column
(*** p<0.001, ** p<0.01, * p<0.05, "" otherwise) for the poster bar fig."""
import sys
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform
from scipy.stats import pearsonr
from gls_tau import tissue_corr_expressed_only, floor_tpm

IN_GTEX_PC_GCT = sys.argv[1]
IN_GTEX_TX_PC_GCT = sys.argv[2]
IN_CHESS_TX_PC_GCT = sys.argv[3]
IN_UBERON_DISTANCE_CSV = sys.argv[4]
OUT_CSV = sys.argv[5]

N_PERM = 9999
RNG = np.random.default_rng(0)
FLOORS = [1.0, 10.0]
METHOD = "pearson"
SOURCE_ORDER = ["gtex_protein_coding", "gtex_v11_tx_protein_coding", "chess_tx_protein_coding"]
SOURCE_PATH = {
    "gtex_protein_coding": IN_GTEX_PC_GCT,
    "gtex_v11_tx_protein_coding": IN_GTEX_TX_PC_GCT,
    "chess_tx_protein_coding": IN_CHESS_TX_PC_GCT,
}
KIND = {"gtex_protein_coding": "gene_centric", "gtex_v11_tx_protein_coding": "tx_centric",
        "chess_tx_protein_coding": "tx_centric"}


def stars(p):
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return ""


uberon = pd.read_csv(IN_UBERON_DISTANCE_CSV, index_col=0)

rows = []
for source in SOURCE_ORDER:
    df = pd.read_csv(SOURCE_PATH[source], sep="\t", skiprows=2)
    tissues = list(df.columns[2:])
    raw_mat = df[tissues].to_numpy(dtype=float)
    for floor in FLOORS:
        mat = floor_tpm(raw_mat, floor=floor)
        R = tissue_corr_expressed_only(mat, method=METHOD)
        dist = 1 - R
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2
        Z = linkage(squareform(dist, checks=False), method="average")
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

        name = f"{source}__floor{int(floor)}"
        rows.append({
            "scenario": name, "source": source, "kind": KIND[source], "floor_threshold": floor,
            "n_common_tissues": n, "mantel_r": r_obs, "p_value": p_value, "stars": stars(p_value),
        })
        print(f"{name}: n_common={n} mantel_r={r_obs:.4f} p={p_value:.4f} {stars(p_value)}", file=sys.stderr)

pd.DataFrame(rows).to_csv(OUT_CSV, index=False)
print(f"-> {OUT_CSV}", file=sys.stderr)
