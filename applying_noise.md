# Noise robustness of tau_unweighted vs tau_weighted

## What we did
Monte Carlo perturbation test: add iid N(0, 0.5) noise to log10(TPM+1) for every
currently-nonzero (TPM>=1, post floor) entry of the full base GCT. Entries already
floored to 0 (TPM<1) are left untouched -- no noise applied. Noisy values that drop
below TPM=1 are not re-floored, but the noisy log10(TPM+1) is clipped at 0 so TPM
never goes negative (`gls_tau.apply_nonzero_log_noise`) -- negligible at low sigma but
a real effect once sigma is large enough for a negative draw to push a small nonzero
entry below 0 in log-space. 200 replicates, one fixed seed. The weight vector
(free_cluster_weights, GTEx onebrain-max) is frozen at its already-fitted value, not
refit per replicate -- it's an aggregate statistic over the whole gene set (~19k
genes) and effectively insensitive to per-cell noise at this scale, so refitting per
replicate would be expensive (differential_evolution) for no real benefit.
tau_unweighted and tau_weighted are recomputed per replicate; the per-gene variance
across the 200 replicates is compared (paired, Wilcoxon signed-rank -- null
hypothesis: both estimators have equal variance under the same perturbation).

Scope: GTEx onebrain-max free_cluster only. (An earlier pass at sigma=0.1, run for
both GTEx and Baron, was retired in favor of this single higher-sigma GTEx-only run.)

Code: `gls_tau.apply_nonzero_log_noise`, `scripts/compute_noise_variance.py`,
`scripts/plot_jackknife_variance_compare.py` (generic column names + a variance-kind
label, reused for both jackknife and noise variance). Snakefile rules
`compute_noise_variance_free_cluster` / `plot_noise_variance_compare_free_cluster`,
wired into `rule all`.

## Result
Null hypothesis rejected (p≈0): tau_weighted (free_cluster) has **higher** per-gene
noise variance than tau_unweighted, not lower.
- GTEx onebrain-max, sigma=0.5: median var 0.00147 (unweighted) -> 0.00162 (weighted);
  frac_above=0.713, frac_below=0.236 (per-gene, weighted vs unweighted variance).

Same direction as the earlier sigma=0.1 pass (frac_above=0.692), slightly stronger at
the higher noise level.

Plausible mechanism: free_cluster weighting concentrates weight onto fewer
"representative" tissues per redundant cluster. That's exactly the right thing to do
against *correlated* redundancy, but under *independent* per-tissue noise it reduces
the effective number of tissues averaged over, so idiosyncratic noise in a
heavily-weighted tissue has an outsized effect -- the reverse of what happens when a
truly redundant tissue is dropped (jackknife variance, tested earlier: weighted was
lower for GTEx). Weighting for redundancy and robustness to independent measurement
noise are not the same objective, and can trade off.

## Next steps for more realistic noise
1. **Heteroscedastic noise.** Real RNA-seq technical noise is Poisson/negative-binomial:
   CV shrinks with expression (~1/sqrt(count) at the Poisson limit), not constant. We
   discussed sigma ~ 1/sqrt(TPM+1) as a heuristic proxy for this -- untested here,
   flat sigma was used for simplicity. Directionally this would push more noise onto
   low/near-floor genes and less onto high-expression genes.
2. **TPM vs counts caveat.** True Poisson noise lives on read counts, not TPM (length-
   and depth-normalized). We don't have raw counts for GTEx's median-TPM-per-tissue
   panel, so any count-based noise model here is a proxy at best.
3. **Correlated (not independent) noise.** Real technical noise often has a batch/library
   component shared across tissues profiled together, which would behave more like the
   jackknife perturbation (correlated) than like this test (fully independent per
   cell) -- could partially explain why this result's direction differs from the
   jackknife result.
4. **Refit weights every K replicates (not every replicate)** to spot-check the
   "weights don't move under noise" assumption cheaply, without paying full DE cost
   every time.
5. **Extend to diag_reg_inv_corr and margin** (frozen weights, same script already
   supports arbitrary WEIGHT_COL/WEIGHTS_CSV) -- not run yet.
6. **Extend to Baron** -- retired when this test moved to sigma=0.5 GTEx-only; could
   be re-added in parallel if wanted.
