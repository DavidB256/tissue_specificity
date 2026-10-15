"""Pairs of transcripts of the same gene that are each (i) isoform_fraction
(tx_total/gene_total TPM, summed across all tissues, tx_total>0 only) >= 1%
and (ii) tx_tau_weighted >= 0.85. For each pair, Wasserstein-1 distance
between the two transcripts' tissue PMFs (per-tissue TPM / tx_total,
same convention as compute_tx_info_content.py's tx_pmf) -- tissues treated
as ordered by their fixed GCT column index (same order for every transcript,
so both PMFs share support positions 0..n_tissues-1)."""
import sys
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import wasserstein_distance

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
OUT_CSV = sys.argv[3]

ISOFORM_FRACTION_THRESHOLD = 0.01
TAU_THRESHOLD = 0.85

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = gct[tissues].to_numpy(dtype=float)
tx_ids = gct["Name"].to_numpy()
gene_ids = gct["Description"].to_numpy()

tx_total = mat.sum(axis=1)
keep = tx_total > 0
mat, tx_ids, gene_ids, tx_total = mat[keep], tx_ids[keep], gene_ids[keep], tx_total[keep]

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(
    columns={"gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted"}
)[["transcript_id", "tx_tau_weighted"]].set_index("transcript_id")["tx_tau_weighted"]

df = pd.DataFrame({"transcript_id": tx_ids, "gene_id": gene_ids, "tx_total": tx_total})
df["gene_total"] = df.groupby("gene_id")["tx_total"].transform("sum")
df["isoform_fraction"] = df["tx_total"] / df["gene_total"]
df["tx_tau_weighted"] = df["transcript_id"].map(tx_tau)
df = df.dropna(subset=["tx_tau_weighted"])

qualifying = df[(df["isoform_fraction"] >= ISOFORM_FRACTION_THRESHOLD) & (df["tx_tau_weighted"] >= TAU_THRESHOLD)]

row_idx = {tx_id: i for i, tx_id in enumerate(tx_ids)}
positions = np.arange(len(tissues), dtype=float)

rows = []
for gene_id, group in qualifying.groupby("gene_id"):
    if len(group) < 2:
        continue
    for (_, a), (_, b) in combinations(group.iterrows(), 2):
        pmf_a = mat[row_idx[a["transcript_id"]]] / a["tx_total"]
        pmf_b = mat[row_idx[b["transcript_id"]]] / b["tx_total"]
        wd = wasserstein_distance(positions, positions, u_weights=pmf_a, v_weights=pmf_b)
        rows.append({
            "gene_id": gene_id,
            "transcript_id_1": a["transcript_id"], "transcript_id_2": b["transcript_id"],
            "isoform_fraction_1": a["isoform_fraction"], "isoform_fraction_2": b["isoform_fraction"],
            "tau_weighted_1": a["tx_tau_weighted"], "tau_weighted_2": b["tx_tau_weighted"],
            "wasserstein_distance": wd,
        })

out = pd.DataFrame(rows).sort_values("wasserstein_distance", ascending=False)
out.to_csv(OUT_CSV, index=False)

print(f"n_qualifying_tx={len(qualifying)} n_genes_with_pairs={out['gene_id'].nunique() if len(out) else 0} "
      f"n_pairs={len(out)} -> {OUT_CSV}", file=sys.stderr)
