# Noise robustness of tau_unweighted vs tau_weighted

## What we did
Monte Carlo perturbation test: add iid N(0, 0.1) noise to log10(TPM+1) for every
currently-nonzero (TPM>=1, post floor) entry of the full base GCT. Entries already
floored to 0 (TPM<1) are left untouched -- no noise applied. Noisy values that drop
below TPM=1 are not re-floored. 200 replicates, one fixed seed. Weight vectors
(free_cluster_weights, max-agg, each source's headline leg) are frozen at their
already-fitted value, not refit per replicate -- each is an aggregate statistic over
the whole gene set (10-20k genes) and effectively insensitive to per-cell noise at
this scale, so refitting per replicate would be expensive (differential_evolution)
for no real benefit. tau_unweighted and tau_weighted are recomputed per replicate; the
per-gene variance across the 200 replicates is compared (paired, Wilcoxon
signed-rank -- null hypothesis: both estimators have equal variance under the same
perturbation).

Ran for both sources at max-agg, headline leg (GTEx onebrain, Baron oneislets).

New code: `gls_tau.apply_nonzero_log_noise`, `scripts/compute_noise_variance.py`,
generalized `scripts/plot_jackknife_variance_compare.py` (now takes column
names + a variance-kind label, reused for both jackknife and noise variance).
Snakefile rules `compute_noise_variance_free_cluster` /
`plot_noise_variance_compare_free_cluster`, wired into `rule all`.

## Result
Null hypothesis rejected in both scenarios (p≈0): tau_weighted (free_cluster) has
**higher** per-gene noise variance than tau_unweighted, not lower.
- GTEx onebrain-max: median var 0.00102 (unweighted) -> 0.00113 (weighted)
- Baron oneislets-max: median var 0.00190 (unweighted) -> 0.00208 (weighted)

Plausible mechanism: free_cluster weighting concentrates weight onto fewer
"representative" tissues per redundant cluster. That's exactly the right thing to do
against *correlated* redundancy, but under *independent* per-tissue noise it reduces
the effective number of tissues averaged over, so idiosyncratic noise in a
heavily-weighted tissue has an outsized effect -- the reverse of what happens when a
truly redundant tissue is dropped (jackknife variance, tested earlier: weighted was
mixed, lower for GTEx, higher for Baron). Weighting for redundancy and robustness to
independent measurement noise are not the same objective, and can trade off.

## Next steps for more realistic noise
1. **Heteroscedastic noise.** Real RNA-seq technical noise is Poisson/negative-binomial:
   CV shrinks with expression (~1/sqrt(count) at the Poisson limit), not constant. We
   discussed sigma ~ 1/sqrt(TPM+1) as a heuristic proxy for this -- untested here,
   flat sigma=0.1 was used for simplicity. Directionally this would push more noise
   onto low/near-floor genes and less onto high-expression genes, likely widening the
   gap this result already shows (weighted tau's few high-weight tissues are usually
   NOT the lowest-expressed ones, so this might shrink or reverse the effect -- worth
   checking).
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
   supports arbitrary WEIGHT_COL/WEIGHTS_CSV) -- not run yet, only free_cluster was
   requested this round.
