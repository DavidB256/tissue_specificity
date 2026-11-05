"""List tx with isoform_frac>=MIN_FRAC and arctanh tau_w diff > MIN_DIFF (multi-tx genes, onebrain weights).
Usage: GCT TX_TAU_CSV WEIGHTS_CSV OUT_CSV [MIN_FRAC=0.1] [MIN_DIFF=1]"""
import sys
import numpy as np
import pandas as pd
from gls_tau import weighted_tau

gct_path, tau_csv, w_csv, out = sys.argv[1:5]
MIN_FRAC = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1
MIN_DIFF = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
SCALE = 0.9999

gct = pd.read_csv(gct_path, sep="\t", skiprows=2)
t = list(gct.columns[2:])
tot = gct[t].sum(axis=1)
k = tot > 0
tx = pd.DataFrame({"transcript_id": gct["Name"][k], "gene_id": gct["Description"][k], "tx_total": tot[k]})
tx["isoform_fraction"] = tx["tx_total"] / tx.groupby("gene_id")["tx_total"].transform("sum")
tau = pd.read_csv(tau_csv).rename(columns={"gene": "transcript_id", "tau_free_weighted": "tx_tau_w"})[
    ["transcript_id", "tx_tau_w"]]
w = pd.read_csv(w_csv).set_index("tissue").loc[t]["weight"].to_numpy()
gm = gct.loc[k, t].groupby(gct.loc[k, "Description"]).sum()
g = pd.DataFrame({"gene_id": gm.index, "gene_tau_w": weighted_tau(gm.to_numpy(), w)})
df = tx.merge(tau, on="transcript_id").merge(g, on="gene_id").dropna(subset=["tx_tau_w", "gene_tau_w"])
df = df[df.groupby("gene_id")["transcript_id"].transform("count") >= 2].copy()
df["arctanh_diff"] = np.arctanh(SCALE * df["tx_tau_w"]) - np.arctanh(SCALE * df["gene_tau_w"])
s = df[(df["isoform_fraction"] >= MIN_FRAC) & (df["arctanh_diff"] > MIN_DIFF)]
s = s.sort_values(["gene_id", "arctanh_diff"], ascending=[True, False])
s.to_csv(out, index=False)
print(f"n_tx={len(s)} n_genes={s['gene_id'].nunique()} -> {out}", file=sys.stderr)
