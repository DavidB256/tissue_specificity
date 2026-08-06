"""Grid-search hard partition of one tissue's genes into 'exp' vs 'beta' groups,
using only that tissue's log10(TPM) and each gene's (tissue-agnostic) tau --
reconstructing the tissue's exp_beta mixture fit from a simple, interpretable rule
instead of the continuous EM assignment.

Decision rule: a linear boundary at an arbitrary angle in the rank-normalized
(log10TPM, tau) plane -- u,v are each gene's percentile rank (0..1) among this
tissue's TPM>=1 genes, so the two axes are on a comparable scale and an "angle"
means something. For angle theta and offset rho: beta group = {u*cos(theta) +
v*sin(theta) >= rho}. Grid over theta (full circle) x rho (quantiles of the
projection) -- looping theta over the full circle at fixed side-convention already
covers both orientations of every physical boundary.

Each grid cell: exp component fit is closed-form MLE (rate = 1/mean, exact). Beta
component fit is a single-start bounded MLE (mixture_core's beta branch, a,b in
[1,50]) via fit_mixture with k=1. Scored by BIC using n_params=6 (1 rate + 4 beta
params + 1 mixing-weight param) -- the same n_params as the continuous exp_beta
model, so the two are directly BIC-comparable."""
import sys
import warnings
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from mixture_core import fit_mixture, logpdf

warnings.filterwarnings("ignore", category=RuntimeWarning)

TISSUE = sys.argv[1]
IN_GCT = sys.argv[2]
IN_TAU = sys.argv[3]
OUT_SUMMARY = sys.argv[4]
OUT_LABELS = sys.argv[5]

N_ANGLES = 12            # 30-degree steps over the full circle
OFFSET_QS = np.arange(0.05, 1.0, 0.05)  # quantiles 0.05, 0.10, .. 0.95 of the projected coordinate
MIN_GROUP_FRAC = 0.02      # both groups must have >=2% of the tissue's genes

gct = pd.read_csv(IN_GCT, sep="\t", skiprows=2, usecols=["Name", TISSUE])
tau_s = pd.read_csv(IN_TAU, usecols=["gene", "tau"]).set_index("gene")["tau"]

tpm = gct[TISSUE].values.astype(float)
valid = tpm >= 1
genes = gct["Name"].values[valid]
x = np.log10(tpm[valid])
tau = tau_s.reindex(genes).values
has_tau = np.isfinite(tau)
genes, x, tau = genes[has_tau], x[has_tau], tau[has_tau]
n = len(x)

u = rankdata(x) / n
v = rankdata(tau) / n

angles = np.linspace(0, 2 * np.pi, N_ANGLES, endpoint=False)

best = None
for theta in angles:
    proj = u * np.cos(theta) + v * np.sin(theta)
    for off in np.quantile(proj, OFFSET_QS):
        is_beta = proj >= off
        n_beta = int(is_beta.sum())
        n_exp = n - n_beta
        if n_beta < MIN_GROUP_FRAC * n or n_exp < MIN_GROUP_FRAC * n:
            continue

        x_exp = x[~is_beta]
        x_beta = x[is_beta]

        rate = 1.0 / x_exp.mean()
        ll_exp = logpdf(x_exp, "expon", {"rate": rate}).sum()

        beta_fit = fit_mixture(x_beta, ["beta"], n_init=1)
        ll_beta = beta_fit["loglik"]

        w_exp, w_beta = n_exp / n, n_beta / n
        total_ll = ll_exp + n_exp * np.log(w_exp) + ll_beta + n_beta * np.log(w_beta)
        bic = 6 * np.log(n) - 2 * total_ll

        if best is None or bic < best["bic"]:
            p = beta_fit["params"][0]
            best = {
                "tissue": TISSUE, "n": n, "theta_deg": float(np.degrees(theta)), "offset": float(off),
                "n_exp": n_exp, "n_beta": n_beta,
                "exp_rate": rate, "beta_a": p["a"], "beta_b": p["b"], "beta_loc": p["loc"], "beta_scale": p["scale"],
                "loglik": total_ll, "bic": bic, "is_beta": is_beta,
            }

labels = pd.DataFrame({"tissue": TISSUE, "gene": genes, "rank_tpm": u, "rank_tau": v,
                        "group": np.where(best["is_beta"], "beta", "exp")})
labels.to_csv(OUT_LABELS, index=False)

summary = {k: v for k, v in best.items() if k != "is_beta"}
pd.DataFrame([summary]).to_csv(OUT_SUMMARY, index=False)
print(f"{TISSUE}: n={n} n_exp={best['n_exp']} n_beta={best['n_beta']} bic={best['bic']:.1f}", file=sys.stderr)
