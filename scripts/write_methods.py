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
lines.append("| model | components | n_params | notes |")
lines.append("|---|---|---|---|")

NOTES = {
    "gauss2": "baseline mixture, no explicit near-zero handling",
    "gauss3": "baseline mixture, no explicit near-zero handling",
    "exp_gauss": "exponential absorbs truncation-boundary pileup; single bulk mode",
    "exp_gauss2": "exponential + 2 gaussians for bulk/tail structure",
    "exp_skewgauss": "skew-normal bulk mode instead of symmetric gaussian",
    "exp_t": "heavy-tailed bulk mode instead of gaussian",
    "exp_gauss_t": "gaussian bulk + separate heavy-tailed component",
    "exp_skewgauss_t": "skew-normal bulk + separate heavy-tailed component",
}

for model, comps in MODELS.items():
    comp_str = " + ".join(COMP_DESC[c] for c in comps)
    lines.append(f"| {model} | {comp_str} | {n_params(comps)} | {NOTES[model]} |")

lines.append(
    "\nFit tables: `results/{model}_fits.csv` (per tissue: aic, bic, ks_stat, ks_pval, "
    "n_params, weights/params as JSON). Cross-model comparison: "
    "`results/model_comparison.csv` (rows=tissue, cols=aic/bic/ks per model). "
    "Grid figures: `figures/{model}_grid.png`, 3x3 panels at BIC quantiles "
    "0,12.5,...,100% across tissues.\n"
)

with open(OUT_MD, "w") as f:
    f.write("\n".join(lines))
