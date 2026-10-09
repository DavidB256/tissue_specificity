# tauQTL detection plan

## Goal
Find genetic variants (cis) assoc. w/ per-donor variation in gene tau
(tissue-specificity). Analogous to cis-eQTL mapping but phenotype = tau, not
expression level.

## Blocker (checked in this session)
No genotype data in `/ccb/salz4-4/dbass13/expr_data/`. GTEx WGS genotypes
are dbGaP/AnVIL controlled-access, not public like the TPM/GCT files already
here. **Cannot proceed past step 1 w/o this.** See unresolved Qs.

## Step 1: per-donor tau phenotype
- Raw per-sample TPM already local: `GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_
  transcripts_tpm.txt.gz` + `sample_tissue_map.csv` (48231 samples).
- Donor ID = sample_id first 2 fields (e.g. `BMS-X4LF`). 1014 donors, median
  22/70 tissues sampled per donor (range 1-41).
- Per donor: build genes x tissues matrix from just their own samples, run
  existing `weighted_tau` (`scripts/gls_tau.py`) w/ hierclust0.1 cluster
  weights -> one tau_weighted per gene per donor.
- **Comparability problem**: donors cover different tissue subsets; tau's
  scale shifts mechanically w/ n_tissues and tissue identity, not just
  biology. Mitigation options (pick one, see unresolved Qs):
  (a) restrict to donors w/ >=K of a fixed "core" tissue panel (most-sampled
  tissues), exclude rest
  (b) include n_tissues_sampled + missing-tissue-pattern PCs as covariates
  in step 3 regression
  (c) both

## Step 2: cis-window + genotype
- Need per-gene TSS -> reuse `extract_transcript_tss.py` output, canonical
  TSS = MANE transcript's TSS (MANE GTF already linked, `link_mane_gtf`
  rule) where available, else widest-used GENCODE transcript.
- Cis window: +/-1Mb around TSS (GTEx eQTL convention).
- Genotype: dosage (0/1/2 or imputed) at all variants (MAF>=0.01, say) in
  window, per donor. **This is the blocked part** -- need GTEx v11 WGS VCF
  (dbGaP) or equivalent array/imputed genotypes, subset to the 1014 donors
  above via donor ID match.

## Step 3: association test
- Per gene: linear regression donor_tau ~ genotype_dosage + covariates,
  one model per cis-variant.
- Covariates: genotype PCs (population stratification), sex, sequencing
  batch/platform (from SampleAttributesDS -- but these are sample-level,
  need donor-level summary or restrict to one reference sample per donor
  per covariate), donor-level tau nuisance covariates from step 1.
- Quantile-normalize tau across donors before regression (tau is bounded
  [0,1], skewed) -- standard eQTL practice, same rationale.
- Tooling: tensorQTL (GPU-friendly, standard GTEx eQTL tool) if genotypes
  arrive as VCF/PLINK; else plain per-gene OLS loop (statsmodels/numpy) if
  variant count small.

## Step 4: multiple testing
- Per gene: permutation (shuffle donor tau labels) to get beta-approximated
  empirical p for best cis-variant (FastQTL/tensorQTL approach).
- Across genes: FDR (Storey q-value or BH) on best-variant p per gene.

## Step 5: validation
- Check top tauQTL genes against known cis-eQTL/sQTL hits (GTEx eQTL
  catalog) for overlap -- a tauQTL overlapping a strong eQTL/sQTL at the
  same locus is more interpretable (e.g. variant shifts isoform usage ->
  shifts tissue specificity) than an isolated hit.
- Sanity check: permute genotype-donor labels genome-wide, confirm null
  p-value distribution is ~uniform (calibration check for the per-donor-tau
  pipeline, independent of real genotype effects).

## Snakemake integration (once genotypes available)
New `scripts/tauqtl_*.py` (plain sys.argv, stderr logging, no `run:` blocks
per CLAUDE.md), new rules under `GENE_SPEC_GTEX_V11_TX_PC`-style naming:
`compute_donor_tau`, `build_gene_cis_windows`, `extract_cis_genotypes`,
`run_tauqtl_mapping`, `tauqtl_permutation_fdr`.

## Unresolved questions
- Genotype access: do you have dbGaP/AnVIL GTEx WGS access? if not, blocked.
- Core-tissue-panel vs covariate approach for donor-coverage issue (1a/1b/1c)?
- min donors per tissue / min tissues per donor cutoffs?
- gene-level or tx-level tau as the QTL phenotype (this pipeline has both)?
- MAF cutoff + cis window size (1Mb std, smaller/larger)?
- quantile-normalize tau, or some other transform?
