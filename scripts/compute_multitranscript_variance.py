"""For each multi-transcript gene (>=2 protein-coding CHESS transcripts in the
gtex_txs panel), compute the range and variance (population, ddof=0) among its
transcripts' tau_unweighted and max_tpm (max TPM across the 9 tissues, no flooring --
every row is a CHESS-annotated transcript, not a candidate call, see
compute_tau_unweighted.py -- the "max agg" scalar used here in place of a per-tissue
breakdown). tau is NaN for transcripts with exactly 0 TPM in every tissue (max<=0);
range/var over tau drop those NaNs and are left NaN themselves if fewer than 2
transcripts remain defined."""
import sys
import numpy as np
import pandas as pd

IN_GCT = sys.argv[1]
IN_TAU_CSV = sys.argv[2]
OUT_CSV = sys.argv[3]

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
mat = gct[tissues].values.astype(float)
max_tpm = mat.max(axis=1)

tau = pd.read_csv(IN_TAU_CSV)[["gene", "tau_unweighted"]]

df = pd.DataFrame({
    "transcript_id": gct["Name"].values,
    "gene_id": gct["Description"].values,
    "max_tpm": max_tpm,
})
df = df.merge(tau, left_on="transcript_id", right_on="gene", how="left").drop(columns="gene")

rows = []
for gene_id, grp in df.groupby("gene_id"):
    n = len(grp)
    if n < 2:
        continue
    tau_valid = grp["tau_unweighted"].dropna().to_numpy()
    tau_range = tau_valid.max() - tau_valid.min() if len(tau_valid) >= 2 else np.nan
    tau_var = tau_valid.var(ddof=0) if len(tau_valid) >= 2 else np.nan
    mt = grp["max_tpm"].to_numpy()
    rows.append({
        "gene_id": gene_id,
        "n_transcripts": n,
        "n_tau_valid": len(tau_valid),
        "tau_range": tau_range,
        "tau_var": tau_var,
        "max_tpm_range": mt.max() - mt.min(),
        "max_tpm_var": mt.var(ddof=0),
    })

out = pd.DataFrame(rows).sort_values("max_tpm_var", ascending=False)
out.to_csv(OUT_CSV, index=False)

print(f"multi-transcript genes={len(out)} (of {df['gene_id'].nunique()} total genes) -> {OUT_CSV}",
      file=sys.stderr)
