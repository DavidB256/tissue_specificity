"""OLS: arctanh(tx_tau_weighted) - arctanh(gene_tau_weighted) ~ 1 +
(tx_length - gene_exonic_length) + gene_tau_weighted + gene_total_tpm +
tx_total_tpm. tau scaled by 0.9999 before arctanh (arctanh(1) undefined).
gene_total_tpm/tx_total_tpm = TPM summed across all tissues. All quantified
transcripts (tx_total>0) of multi-tx genes (n_tx>=2), not just
"significant" ones."""
import sys
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

IN_GCT = sys.argv[1]
IN_TX_TAU_CSV = sys.argv[2]
IN_INFO_CONTENT_CSV = sys.argv[3]
IN_LENGTH_FRACTION_CSV = sys.argv[4]
OUT_CSV = sys.argv[5]

SCALE = 0.9999


def arctanh_tau(tau):
    return np.arctanh(SCALE * tau)


gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(gct.columns[2:])
tx_total = gct[tissues].sum(axis=1)
keep = tx_total > 0

tx = pd.DataFrame({
    "transcript_id": gct["Name"][keep],
    "gene_id": gct["Description"][keep],
    "tx_total_tpm": tx_total[keep],
})
tx["gene_total_tpm"] = tx.groupby("gene_id")["tx_total_tpm"].transform("sum")

tx_tau = pd.read_csv(IN_TX_TAU_CSV).rename(columns={
    "gene": "transcript_id", "tau_free_weighted": "tx_tau_weighted",
})[["transcript_id", "tx_tau_weighted"]]
info = pd.read_csv(IN_INFO_CONTENT_CSV)[["gene_id", "gene_tau_weighted"]]
length = pd.read_csv(IN_LENGTH_FRACTION_CSV)[["transcript_id", "tx_length", "gene_exonic_length"]]

df = tx.merge(tx_tau, on="transcript_id", how="inner")
df = df.merge(info, on="gene_id", how="inner")
df = df.merge(length, on="transcript_id", how="inner")
df = df.dropna(subset=["tx_tau_weighted", "gene_tau_weighted", "tx_length", "gene_exonic_length",
                        "gene_total_tpm", "tx_total_tpm"])

n_tx_per_gene = df.groupby("gene_id")["transcript_id"].transform("count")
df = df[n_tx_per_gene >= 2].copy()

df["deviation_arctanh"] = arctanh_tau(df["tx_tau_weighted"]) - arctanh_tau(df["gene_tau_weighted"])
df["length_diff"] = df["tx_length"] - df["gene_exonic_length"]

model = smf.ols(
    "deviation_arctanh ~ 1 + length_diff + gene_tau_weighted + gene_total_tpm + tx_total_tpm",
    data=df,
).fit()

summary_text = model.summary().as_text()
print(summary_text, file=sys.stderr)

out = pd.DataFrame({
    "term": model.params.index,
    "coef": model.params.values,
    "std_err": model.bse.values,
    "t": model.tvalues.values,
    "p_value": model.pvalues.values,
})
out.attrs["r_squared"] = model.rsquared
out.to_csv(OUT_CSV, index=False)

print(f"\nn_tx={len(df)} r_squared={model.rsquared:.4f} adj_r_squared={model.rsquared_adj:.4f} -> {OUT_CSV}",
      file=sys.stderr)
