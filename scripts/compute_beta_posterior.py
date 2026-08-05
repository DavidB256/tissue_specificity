"""Per gene-tissue posterior probability of belonging to the exp_beta model's beta
(vs exponential) component: P(beta|x) = w2*beta_pdf(x) / (w1*expon_pdf(x) + w2*beta_pdf(x)),
evaluated at each gene's log10(TPM) -- the domain the model was fit on (log10(TPM>=1)
per tissue). TPM<1 falls outside that domain -> NaN. Adds beta_posterior_mean and
beta_posterior_median across tissues (bounded in [0,1], so unlike exp_beta_ratio no
outlier blow-up -- mean is well-behaved here, median kept alongside for comparison)."""
import sys
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from mixture_core import logpdf

IN_GCT = sys.argv[1]
IN_FITS = sys.argv[2]
OUT_CSV = sys.argv[3]

df = pd.read_csv(IN_GCT, sep="\t", skiprows=2)
tissues = list(df.columns[2:])
genes = df["Name"].values
mat = df[tissues].values.astype(float)

fits = pd.read_csv(IN_FITS).set_index("tissue")

post_cols = {}
for j, t in enumerate(tissues):
    row = fits.loc[t]
    x = mat[:, j]
    valid = x >= 1
    with np.errstate(divide="ignore", invalid="ignore"):
        log10x = np.log10(x)
        expon_lp = np.log(row["weight_1"]) + logpdf(log10x, "expon", {"rate": row["comp1_rate"]})
        beta_lp = np.log(row["weight_2"]) + logpdf(log10x, "beta", {"a": row["comp2_a"], "b": row["comp2_b"],
                                                                       "loc": row["comp2_loc"], "scale": row["comp2_scale"]})
        denom = logsumexp(np.stack([expon_lp, beta_lp]), axis=0)
        post = np.exp(beta_lp - denom)
    post = np.where(valid, post, np.nan)
    post_cols[t] = post

out = pd.DataFrame(post_cols)
out.insert(0, "gene", genes)
out["beta_posterior_mean"] = out[tissues].mean(axis=1, skipna=True)
out["beta_posterior_median"] = out[tissues].median(axis=1, skipna=True)
out.to_csv(OUT_CSV, index=False)
print(f"genes={len(out)} tissues={len(tissues)}", file=sys.stderr)
