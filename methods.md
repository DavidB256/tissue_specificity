# Methods: log10(TPM>=1) mixture models

Per tissue, TPM values <1 are dropped, remaining values log10-transformed. Each model below is a finite mixture fit independently per tissue via direct maximum likelihood (`scipy.optimize.minimize`, L-BFGS-B, 5 random restarts, best log-likelihood kept) rather than EM, so skew-normal and Student-t components (no closed-form EM M-step) can be mixed freely with the others. See `scripts/mixture_core.py`.

Component weights are free except for the sum-to-1 constraint (k components -> k-1 free weight parameters). Params/component: norm=2, expon=1, skewnorm=3, t=3, beta=4, gamma=3.

Student-t degrees of freedom is capped at df<=50 (hardcoded clip, not a smooth reparameterization): beyond df~50 a t-distribution is statistically indistinguishable from a Gaussian (excess kurtosis 6/(df-4) < 0.13), so the extra shape parameter buys nothing and an unconstrained fit tends to run away toward df->infinity wherever the data don't need real heavy tails.

Skew-normal shape (a) is capped at |a|<=100, and beta/gamma shape params (a, b) at <=100, for the same reason: an unconstrained fit can run these away to extreme values chasing a degenerate near-Gaussian shape.

| model | components | n_params | notes |
|---|---|---|---|
| gauss2 | Gaussian (mean, std) + Gaussian (mean, std) | 5 | baseline mixture, no explicit near-zero handling |
| exp_gauss | Exponential (rate) + Gaussian (mean, std) | 4 | exponential absorbs truncation-boundary pileup; single bulk mode |
| exp_gauss2 | Exponential (rate) + Gaussian (mean, std) + Gaussian (mean, std) | 7 | exponential + 2 gaussians for bulk/tail structure |
| exp_skewgauss | Exponential (rate) + Skew-normal (loc, scale, shape) | 5 | skew-normal bulk mode instead of symmetric gaussian |
| exp_skewgauss2 | Exponential (rate) + Skew-normal (loc, scale, shape) + Skew-normal (loc, scale, shape) | 9 | 2 skew-normals for bulk/tail structure, no symmetric gaussian |
| exp_skewgauss_gauss | Exponential (rate) + Skew-normal (loc, scale, shape) + Gaussian (mean, std) | 8 | skew-normal bulk + separate symmetric gaussian component |
| exp_t | Exponential (rate) + Student-t (loc, scale, df) | 5 | heavy-tailed (df<=50) bulk mode instead of gaussian |
| exp_gauss_t | Exponential (rate) + Gaussian (mean, std) + Student-t (loc, scale, df) | 8 | gaussian bulk + separate heavy-tailed (df<=50) component |
| exp_skewgauss_t | Exponential (rate) + Skew-normal (loc, scale, shape) + Student-t (loc, scale, df) | 9 | skew-normal bulk + separate heavy-tailed (df<=50) component |
| exp_beta | Exponential (rate) + Beta type I (a, b, loc, scale) | 6 | bounded (4-param) beta bulk mode instead of gaussian |
| exp_gamma | Exponential (rate) + Gamma (a, loc, scale) | 5 | half-bounded gamma bulk mode instead of gaussian |

Fit tables: `model_fits/{model}_fits.csv` (per tissue: aic, bic, ks_stat, ks_pval, n_params, weights/params as JSON). Cross-model comparison: `model_fits/model_comparison.csv` (rows=tissue, cols=aic/bic/ks per model). Grid figures: `distribution_fits/{model}_grid.png`, 3x3 panels at BIC quantiles 0,12.5,...,100% across tissues. Tissue clustering: `tissue_pca/{model}_pca.png` (PC1-vs-PC2, PC1-vs-PC3 scatter of tissues in the model's fitted, standardized parameter space; every tissue labeled).
