"""Per gene-tissue likelihood ratio of the exp_beta model's exponential component
vs its beta component, evaluated at each gene's log10(TPM) -- the domain the model
was fit on (log10(TPM>=1) per tissue). TPM<1 falls outside that domain -> NaN.
Adds exp_beta_ratio = mean across tissues (ignoring NaN), and exp_beta_ratio_median
= median across tissues (more robust to the occasional extreme per-tissue ratio when
a gene's expression sits right at the edge of a tissue's fitted beta support)."""
import sys
import pandas as pd
import numpy as np
from mixture_core import logpdf

IN_GCT = sys.argv[1]
IN_FITS = sys.argv[2]
OUT_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)

fits = pd.read_csv(IN_FITS).set_index("tissue")

ratio_cols = {}
for j, t in enumerate(tissues):
    row = fits.loc[t]
    x = mat[:, j]
    valid = x >= 1
    with np.errstate(divide="ignore", invalid="ignore"):
        log10x = np.log10(x)
        expon_lp = logpdf(log10x, "expon", {"rate": row["comp1_rate"]})
        beta_lp = logpdf(log10x, "beta", {"a": row["comp2_a"], "b": row["comp2_b"],
                                           "loc": row["comp2_loc"], "scale": row["comp2_scale"]})
        ratio = np.exp(expon_lp - beta_lp)
    ratio = np.where(valid & np.isfinite(ratio), ratio, np.nan)
    ratio_cols[t] = ratio

out = pd.DataFrame(ratio_cols)
out.insert(0, "gene", genes)
out["exp_beta_ratio"] = out[tissues].mean(axis=1, skipna=True)
out["exp_beta_ratio_median"] = out[tissues].median(axis=1, skipna=True)
out.to_csv(OUT_CSV, index=False)
print(f"genes={len(out)} tissues={len(tissues)}", file=sys.stderr)
