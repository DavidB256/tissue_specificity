"""Summarize the top-N genes by weighted_kl (TPM-share-weighted average of
KL(tx tissue pmf || gene tissue pmf), compute_tx_info_content.py) for
gtex_v11_tx_protein_coding: gene identity, isoform-count/entropy stats, tau
columns, and -- from the per-gene pickle -- the dominant (highest-TPM) and
most-divergent (highest individual KL) transcript, which are not always the
same isoform. Not wired into snakefile -- ad hoc diagnostic, run manually."""
import sys
import pickle
import pandas as pd

IN_CSV = sys.argv[1]
IN_PICKLE = sys.argv[2]
IN_NAME_MAP_CSV = sys.argv[3]
OUT_CSV = sys.argv[4]
TOP_N = int(sys.argv[5]) if len(sys.argv) > 5 else 20

df = pd.read_csv(IN_CSV)
name_map = pd.read_csv(IN_NAME_MAP_CSV)
name_map["gene_base"] = name_map["gene_id"].str.split(".").str[0]
df["gene_base"] = df["gene_id"].str.split(".").str[0]
df = df.merge(name_map[["gene_base", "gene_name"]], on="gene_base", how="left").drop(columns="gene_base")

with open(IN_PICKLE, "rb") as f:
    kl_dict = pickle.load(f)

top = df.nlargest(TOP_N, "weighted_kl").reset_index(drop=True)

rows = []
for _, row in top.iterrows():
    tx_entries = kl_dict[row["gene_id"]]
    dominant_tx = max(tx_entries, key=lambda t: tx_entries[t][0])
    top_kl_tx = max(tx_entries, key=lambda t: tx_entries[t][1])
    total_tpm = sum(v[0] for v in tx_entries.values())
    rows.append({
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "n_tx": row["n_tx"],
        "H": row["H"],
        "H_norm": row["H_norm"],
        "effective_n": row["effective_n"],
        "weighted_kl": row["weighted_kl"],
        "tx_tau_weighted_mean": row["tx_tau_weighted_mean"],
        "tx_tau_weighted_rms": row["tx_tau_weighted_rms"],
        "gene_tau_weighted": row["gene_tau_weighted"],
        "dominant_tx_id": dominant_tx,
        "dominant_tx_tpm_share": tx_entries[dominant_tx][0] / total_tpm,
        "dominant_tx_kl": tx_entries[dominant_tx][1],
        "top_kl_tx_id": top_kl_tx,
        "top_kl_tx_tpm_share": tx_entries[top_kl_tx][0] / total_tpm,
        "top_kl_tx_kl": tx_entries[top_kl_tx][1],
        "dominant_is_top_kl": dominant_tx == top_kl_tx,
    })

out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)

n_dominant_is_top = int(out["dominant_is_top_kl"].sum())
print(f"top {TOP_N} genes by weighted_kl: range [{out['weighted_kl'].min():.4f}, {out['weighted_kl'].max():.4f}]",
      file=sys.stderr)
print(f"n_tx range: [{out['n_tx'].min()}, {out['n_tx'].max()}], median {out['n_tx'].median():.0f}",
      file=sys.stderr)
print(f"dominant tx is also the highest-individual-KL tx for {n_dominant_is_top}/{TOP_N} genes",
      file=sys.stderr)
print(f"gene_tau_weighted range: [{out['gene_tau_weighted'].min():.4f}, {out['gene_tau_weighted'].max():.4f}]",
      file=sys.stderr)
print(f"-> {OUT_CSV}", file=sys.stderr)
print(file=sys.stderr)
print(out[["gene_name", "n_tx", "weighted_kl", "gene_tau_weighted", "dominant_tx_tpm_share",
           "dominant_is_top_kl"]].to_string(index=False), file=sys.stderr)
