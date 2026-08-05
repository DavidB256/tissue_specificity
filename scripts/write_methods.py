"""Generate methods.md: table of all mixture models, their components, and param counts."""
import sys
from mixture_core import n_params, COMP_NPARAMS
from model_registry import MODELS

OUT_MD = sys.argv[1]

COMP_DESC = {
    "expon": "Exponential (rate)",
    "norm": "Gaussian (mean, std)",
    "skewnorm": "Skew-normal (loc, scale, shape)",
    "t": "Student-t (loc, scale, df)",
    "beta": "Beta type I (a, b, loc, scale)",
    "gamma": "Gamma (a, loc, scale)",
}

lines = []
lines.append("# Methods: log10(TPM>=1) mixture models\n")
lines.append(
    "Per tissue, TPM values <1 are dropped, remaining values log10-transformed. "
    "Each model below is a finite mixture fit independently per tissue via direct "
    "maximum likelihood (`scipy.optimize.minimize`, L-BFGS-B, 5 random restarts, "
    "best log-likelihood kept) rather than EM, so skew-normal and Student-t "
    "components (no closed-form EM M-step) can be mixed freely with the others. "
    "See `scripts/mixture_core.py`.\n"
)
lines.append(
    "Component weights are free except for the sum-to-1 constraint (k components -> "
    "k-1 free weight parameters). Params/component: " +
    ", ".join(f"{k}={v}" for k, v in COMP_NPARAMS.items()) + ".\n"
)
lines.append(
    "Student-t degrees of freedom is capped at df<=50 (hardcoded clip, not a smooth "
    "reparameterization): beyond df~50 a t-distribution is statistically "
    "indistinguishable from a Gaussian (excess kurtosis 6/(df-4) < 0.13), so the extra "
    "shape parameter buys nothing and an unconstrained fit tends to run away toward "
    "df->infinity wherever the data don't need real heavy tails.\n"
)
lines.append(
    "Skew-normal shape (a) is capped at |a|<=100, and gamma shape (a) at <=100, for "
    "the same reason: an unconstrained fit can run these away to extreme values "
    "chasing a degenerate near-Gaussian shape. Beta shape params (a, b) are bounded "
    "to [1, 50]: a,b>=1 keeps the beta concave-down (matches the intended bulk-mode "
    "role), and <=50 bounds the search space / prevents runaway.\n"
)
lines.append("| model | components | n_params | notes |")
lines.append("|---|---|---|---|")

NOTES = {
    "gauss2": "baseline mixture, no explicit near-zero handling",
    "exp_gauss": "exponential absorbs truncation-boundary pileup; single bulk mode",
    "exp_gauss2": "exponential + 2 gaussians for bulk/tail structure",
    "exp_skewgauss": "skew-normal bulk mode instead of symmetric gaussian",
    "exp_skewgauss2": "2 skew-normals for bulk/tail structure, no symmetric gaussian",
    "exp_skewgauss_gauss": "skew-normal bulk + separate symmetric gaussian component",
    "exp_t": "heavy-tailed (df<=50) bulk mode instead of gaussian",
    "exp_gauss_t": "gaussian bulk + separate heavy-tailed (df<=50) component",
    "exp_skewgauss_t": "skew-normal bulk + separate heavy-tailed (df<=50) component",
    "exp_beta": "bounded (4-param) beta bulk mode instead of gaussian",
    "exp_gamma": "half-bounded gamma bulk mode instead of gaussian",
}

for model, comps in MODELS.items():
    comp_str = " + ".join(COMP_DESC[c] for c in comps)
    lines.append(f"| {model} | {comp_str} | {n_params(comps)} | {NOTES[model]} |")

lines.append(
    "\nFit tables: `model_fits/{model}_fits.csv` (per tissue: aic, bic, ks_stat, "
    "ks_pval, n_params, weights/params as JSON). Cross-model comparison: "
    "`model_fits/model_comparison.csv` (rows=tissue, cols=aic/bic/ks per model). "
    "Grid figures: `distribution_fits/{model}_grid.png`, 3x3 panels at BIC quantiles "
    "0,12.5,...,100% across tissues. Tissue clustering: `tissue_pca/{model}_pca.png` "
    "(PC1-vs-PC2, PC1-vs-PC3 scatter of tissues in the model's fitted, standardized "
    "parameter space; every tissue labeled).\n"
)

with open(OUT_MD, "w") as f:
    f.write("\n".join(lines))
