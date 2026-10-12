import sys
sys.path.insert(0, "scripts")
import pandas as pd
from model_registry import MODELS

DATA = "/ccb/salz4-4/dbass13/expr_data"
GCT = f"{DATA}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
MODEL_NAMES = list(MODELS.keys())
TISSUES = list(pd.read_csv(GCT, sep="\t", skiprows=2, nrows=0).columns[2:])

PROCESSED = f"{DATA}/processed"
MODEL_FITS = f"{DATA}/model_fits"
DIST_FIGS = f"{DATA}/distribution_fits"
PCA_FIGS = f"{DATA}/tissue_pca"
GENE_SPEC = f"{DATA}/gene_specificity"

MANE_GTF_SRC = "/ccb/salz4-4/dbass13/recursive_introns_data/annotations/MANE.v1.5.gtf"
GENCODE_V47_GTF = "/ccb/salz4-4/dbass13/recursive_introns_data/annotations/gencode.v47.annotation.gtf"
GENCODE_V47_GTF_LOCAL = f"{DATA}/gencode.v47.annotation.gtf"

# Baron et al. (GSE84133) pancreas scRNA-seq -> per-cell-type TPM matrix. Gene
# length for TPM = MANE_Select summed exon length, falling back to the longest
# GENCODE transcript (by summed exon bp) for genes not in MANE.
BARON_DIR = f"{DATA}/Baron"
BARON_HUMAN_SAMPLES = [
    f"{BARON_DIR}/GSM2230757_human1_umifm_counts.csv.gz",
    f"{BARON_DIR}/GSM2230758_human2_umifm_counts.csv.gz",
    f"{BARON_DIR}/GSM2230759_human3_umifm_counts.csv.gz",
    f"{BARON_DIR}/GSM2230760_human4_umifm_counts.csv.gz",
]
GENE_SPEC_BARON = f"{GENE_SPEC}/baron"
BARON_ONEISLETS_MEMBERS = ["alpha", "beta", "gamma", "delta", "epsilon"]

wildcard_constraints:
    tissue = r"[^/]+"

# ============================================================================
# Protein-coding-only tissue/cell-type-specificity pipeline (GTEx + Baron).
#
# Only protein-coding genes are analyzed going forward (all-genes GTEx/Baron
# analysis retired -- see scripts/old/ and git history for what used to build
# GENE_SPEC/gtex and GENE_SPEC/baron; the underlying raw data is untouched on
# disk, just no longer targeted by this snakefile).
#
# Five decision axes, coordinated through wildcards instead of copy-pasted
# rule blocks:
#   source (gtex_protein_coding | baron_protein_coding) -- SOURCES registry
#   resolution/"leg" (onebrain|oneislets | hierclust0.1) -- per-source, see
#       SOURCES[source]["legs"]
#   combination method (mean | max) -- {agg_suffix} wildcard, "_mean" or
#       "_max", always explicit (no bare/unsuffixed default -- mean and max
#       are on equal footing)
#   gene filter -- fixed to protein_coding (the whole point of this section)
#   weighting method (free_weights | free_cluster_weights | diag_reg_inv_corr |
#       margin) -- free_weights/free_cluster_weights share one optimize/combine
#       script (scripts/optimize_free_weights.py, MODE=tissue|cluster);
#       diag_reg_inv_corr ("diagonal-regularized inverse correlation") is
#       closed-form (compute_weighted_tau.py) from the tissue x tissue
#       correlation matrix; margin (compute_margin_weighted_tau.py) is
#       closed-form directly from the expression matrix, no correlation matrix
#       at all (gls_tau.margin_weights). Both diag_reg_inv_corr and margin are
#       agg/leg-independent (no leg/agg-specific weight fit needed).
# ============================================================================

GENE_SPEC_GTEX_PC = f"{GENE_SPEC}/gtex_protein_coding"
GENE_SPEC_BARON_PC = f"{GENE_SPEC}/baron_protein_coding"

SOURCES = {
    "gtex_protein_coding": dict(
        legs={"onebrain": ("prefix", "Brain"), "hierclust0.1": ("cluster", None)},
        n_free=60, n_free_cluster=20, display="GTEx (protein-coding)",
    ),
    "baron_protein_coding": dict(
        legs={"oneislets": ("group", BARON_ONEISLETS_MEMBERS), "hierclust0.1": ("cluster", None)},
        n_free=10, n_free_cluster=10, display="Baron pancreas (protein-coding)",
    ),
}
STUB = {"gtex_protein_coding": "gene_median_tpm", "baron_protein_coding": "celltype_tpm",
        "gtex_v11_tx": "transcript_median_tpm", "gtex_v11_tx_gene": "gene_median_tpm",
        "gtex_v11_tx_protein_coding": "transcript_median_tpm"}
RAW_GCT = {"gtex_protein_coding": GCT, "baron_protein_coding": f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz"}
# legacy paths (predate this wildcard scheme, kept as literals rather than migrated --
# no data movement, see streamline-snakefile.md plan): GTEx nests clusters under a
# "full/" subdir, Baron doesn't.
CLUSTERS_CSV = {
    "gtex_protein_coding": f"{GENE_SPEC_GTEX_PC}/cluster_weights/full/tissue_clusters.csv",
    "baron_protein_coding": f"{GENE_SPEC_BARON_PC}/cluster_weights/tissue_clusters_hierclust0.1.csv",
    "gtex_v11_tx": f"{GENE_SPEC}/gtex_v11_tx/cluster_weights/tissue_clusters_hierclust0.1.csv",
    "gtex_v11_tx_gene": f"{GENE_SPEC}/gtex_v11_tx_gene/cluster_weights/tissue_clusters_hierclust0.1.csv",
    "gtex_v11_tx_protein_coding": f"{GENE_SPEC}/gtex_v11_tx_protein_coding/cluster_weights/tissue_clusters_hierclust0.1.csv",
}
CLUSTERS_TAU_CSV = {
    "gtex_protein_coding": f"{GENE_SPEC_GTEX_PC}/cluster_weights/full/weighted_tau.csv",
    "baron_protein_coding": f"{GENE_SPEC_BARON_PC}/cluster_weights/weighted_tau_hierclust0.1.csv",
    "gtex_v11_tx": f"{GENE_SPEC}/gtex_v11_tx/cluster_weights/weighted_tau_hierclust0.1.csv",
    "gtex_v11_tx_gene": f"{GENE_SPEC}/gtex_v11_tx_gene/cluster_weights/weighted_tau_hierclust0.1.csv",
    "gtex_v11_tx_protein_coding": f"{GENE_SPEC}/gtex_v11_tx_protein_coding/cluster_weights/weighted_tau_hierclust0.1.csv",
}
TAU_PLOT_TITLE = {
    "gtex_protein_coding": "log10 TPM, TPM<1->0, GTEx tissues (protein-coding)",
    "baron_protein_coding": "log10 TPM, TPM<1->0, Baron pancreas cell types (protein-coding)",
    "gtex_v11_tx": "log10 TPM, TPM<1->0, GTEx v11 transcripts (protein-coding)",
    "gtex_v11_tx_protein_coding": "log10 TPM, TPM<1->0, GTEx v11 transcripts (protein-coding, strict transcript biotype)",
}
# headline free_cluster_weights_max leg per source, for the crossing_thresholds gene lists
CROSSING_THRESHOLD_LEG = {"gtex_protein_coding": "onebrain", "baron_protein_coding": "oneislets"}
NOISE_SIGMA = 0.5
NOISE_N_REPLICATES = 200
NOISE_SEED = 0
CROSSING_CATEGORIES = ["above0.5_to_below", "below0.5_to_above", "above0.8_to_below", "below0.8_to_above"]
NEURAL_CROSSING_THRESHOLD = 0.85


def base_gct(source):
    return f"{GENE_SPEC}/{source}/{STUB[source]}_protein_coding.gct.gz"


# ---- gtex_v11_tx: RSEM transcript-level tau on the full 68-tissue GTEx v11 panel
# (Ensembl transcript/gene IDs). Downstream analysis (tau/multi-transcript-variance/
# figures) is wildcarded over {txs_source} below (TXS_BASE_GCT/TXS_META_CSV) --
# currently a single source, kept wildcarded for easy extension to a future
# transcript-level source. TPM<1 is floored immediately after the per-tissue median,
# in build_gtex_v11_tx_gct.py. ----

GENE_SPEC_GTEX_V11_TX = f"{GENE_SPEC}/gtex_v11_tx"
GTEX_V11_TX_TPM_URL = "https://storage.googleapis.com/adult-gtex/bulk-gex/v11/rna-seq/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz"
GTEX_V11_TX_TPM_RAW = f"{DATA}/GTEx_Analysis_2025-08-22_v11_RSEMv1.3.3_transcripts_tpm.txt.gz"
GTEX_V11_SAMPLE_ATTRS_URL = "https://storage.googleapis.com/adult-gtex/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
GTEX_V11_SAMPLE_ATTRS_RAW = f"{DATA}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt"
GTEX_V11_TX_BASE_GCT = f"{GENE_SPEC_GTEX_V11_TX}/transcript_median_tpm_protein_coding.gct.gz"

# ---- gtex_v11_tx_protein_coding: gtex_v11_tx further filtered to transcripts whose
# OWN GENCODE v47 transcript_type=="protein_coding" (stricter than gtex_v11_tx's
# gene-level protein_coding filter, which still admits e.g. NMD/retained-intron
# isoforms of a coding gene). ----

GENE_SPEC_GTEX_V11_TX_PC = f"{GENE_SPEC}/gtex_v11_tx_protein_coding"
GTEX_V11_TX_PC_BASE_GCT = f"{GENE_SPEC_GTEX_V11_TX_PC}/transcript_median_tpm_protein_coding.gct.gz"

TXS_SOURCES_RE = r"gtex_v11_tx|gtex_v11_tx_protein_coding"
TXS_BASE_GCT = {"gtex_v11_tx": GTEX_V11_TX_BASE_GCT, "gtex_v11_tx_protein_coding": GTEX_V11_TX_PC_BASE_GCT}
TXS_META_CSV = {
    "gtex_v11_tx": f"{GENE_SPEC_GTEX_V11_TX}/transcript_meta.csv",
    "gtex_v11_tx_protein_coding": f"{GENE_SPEC_GTEX_V11_TX_PC}/transcript_meta.csv",
}

# ---- gtex_v11_tx DE free_cluster_weights, tx and gene resolution, hierclust0.1 leg
# only, max-agg only (per explicit choice -- see commit history). Gene-level GCT
# (gtex_v11_tx_gene) sums the already-floored transcript GCT per gene_id; reuses
# base_gct()/STUB like any other source. Tissue clustering (hierclust0.1) is computed
# once from the tx-level correlation (already built for the corr_and_inverse figure)
# and shared by both resolutions -- clustering is a tissue-tissue relationship,
# independent of whether the expression matrix is at tx or gene resolution. Reuses the
# generic build_cluster_leg_gct/compute_weighted_tau_leg/optimize_free_cluster_weights_chain
# rules (source-wildcarded, no SOURCES-dict coupling); combine_free_cluster_weights_chains
# is NOT reused since it reads SOURCES[source]["n_free_cluster"] -- a bespoke combine
# rule below hardcodes N_FREE_CLUSTER_GTEX_V11_TX instead, to avoid adding gtex_v11_tx*
# to SOURCES (which would also auto-wire unrequested margin/diag_reg_inv_corr/
# noise-variance/free_weights targets via rule all's source=SOURCES expansions). ----

GTEX_V11_TX_GENE_BASE_GCT = base_gct("gtex_v11_tx_gene")
GTEX_V11_TX_RES_RE = r"gtex_v11_tx|gtex_v11_tx_gene|gtex_v11_tx_protein_coding"
N_FREE_CLUSTER_GTEX_V11_TX = 20
# which dir's own tissue_correlation.csv/linkage.csv each DE "res" clusters against --
# gtex_v11_tx_gene reuses the tx-level basis (same tissue-tissue relationship, see
# comment above); gtex_v11_tx_protein_coding gets its own (different transcript set).
CLUSTER_CORR_BASIS = {
    "gtex_v11_tx": "gtex_v11_tx",
    "gtex_v11_tx_gene": "gtex_v11_tx",
    "gtex_v11_tx_protein_coding": "gtex_v11_tx_protein_coding",
}


def leg_gct_path(source, leg, agg_suffix=""):
    return f"{GENE_SPEC}/{source}/{leg}{agg_suffix}/{STUB[source]}_{leg}.gct.gz"


def leg_tau_path(source, leg, agg_suffix=""):
    return f"{GENE_SPEC}/{source}/{leg}{agg_suffix}/weighted_tau.csv"


def leg_members(source, leg):
    kind, val = SOURCES[source]["legs"][leg]
    if kind == "prefix":
        return [t for t in TISSUES if t.split("_")[0] == val]
    return val  # kind == "group": explicit member list


LEGS_RE = r"onebrain|oneislets|hierclust0\.1"
AGGS = ["mean", "max"]
AGG_SUFFIX_RE = r"_mean|_max"

wildcard_constraints:
    source = r"gtex_protein_coding|baron_protein_coding|gtex_v11_tx|gtex_v11_tx_gene|gtex_v11_tx_protein_coding",
    leg = LEGS_RE,
    agg_suffix = AGG_SUFFIX_RE,
    method = r"free_weights|free_cluster_weights"

# figure-path convention is uniform across source/method/agg (no more legacy
# flat-dir special case -- mean and max nest under {method}{agg_suffix}/figures/
# identically, see plot_free_weights_baron_legacy removal in git history).


def figures_dir(source, method, agg_suffix):
    return f"{GENE_SPEC}/{source}/{method}{agg_suffix}/figures"


def all_free_weights_targets(kind):
    """kind: 'weights' | 'chains_summary' -- combine_*_chains outputs, for every
    (source, method, leg, agg) combo (backfilled to full coverage, see
    streamline-snakefile.md plan's 'backfill' decision)."""
    return [f"{GENE_SPEC}/{source}/{method}_{agg}/{leg}/{kind}.csv"
            for source, cfg in SOURCES.items()
            for method in ("free_weights", "free_cluster_weights")
            for leg in cfg["legs"]
            for agg in AGGS]


rule all:
    input:
        f"{MODEL_FITS}/model_comparison.csv",
        "methods.md",
        expand(f"{DIST_FIGS}/{{model}}_grid.png", model=MODEL_NAMES),
        expand(f"{PCA_FIGS}/{{model}}_pca.png", model=MODEL_NAMES),
        f"{MODEL_FITS}/pca_compressibility_comparison.csv",
        f"{MODEL_FITS}/pca_compressibility_comparison_allgenes.csv",
        # per-source: correlation heatmap, dendrogram, full-resolution cluster weights
        expand(f"{GENE_SPEC}/{{source}}/figures/tissue_correlation_clustermap.png", source=SOURCES),
        expand(f"{GENE_SPEC}/{{source}}/figures/corr_and_inverse.png", source=SOURCES),
        expand(f"{GENE_SPEC}/{{source}}/figures/tissue_clusters_dendrogram_hierclust0.1.png", source=SOURCES),
        # per-source, per-leg, per-agg: leg weighted_tau (reference tau for free_weights)
        [leg_tau_path(source, leg, f"_{agg}")
         for source, cfg in SOURCES.items() for leg in cfg["legs"] for agg in AGGS],
        # free_weights / free_cluster_weights: chains_summary + weights + tau, all combos
        all_free_weights_targets("weights"),
        all_free_weights_targets("chains_summary"),
        # weights bar plots, all combos (mean and max on equal footing)
        [figures_dir(source, method, f"_{agg}") + f"/tissue_weight_free{'_cluster' if method == 'free_cluster_weights' else ''}_{leg}.png"
         for source, cfg in SOURCES.items() for method in ("free_weights", "free_cluster_weights")
         for leg in cfg["legs"] for agg in AGGS],
        # tau-vs-reference comparison plots, all combos
        [figures_dir("gtex_protein_coding", method, f"_{agg}") +
         (f"/tau_free_vs_full_{leg}.png" if method == "free_weights" else f"/tau_free_cluster_vs_full_{leg}.png")
         for method in ("free_weights", "free_cluster_weights")
         for leg in SOURCES["gtex_protein_coding"]["legs"] for agg in AGGS],
        [figures_dir("baron_protein_coding", method, f"_{agg}") +
         (f"/tau_scenarios_scatter_{leg}.png" if method == "free_weights" else f"/tau_free_cluster_scenarios_scatter_{leg}.png")
         for method in ("free_weights", "free_cluster_weights")
         for leg in SOURCES["baron_protein_coding"]["legs"] for agg in AGGS],
        # genes crossing tau 0.5/0.8 thresholds, headline free_cluster_weights_max leg per source
        [f"{GENE_SPEC}/{source}/free_cluster_weights_max/crossing_thresholds/{leg}_{cat}.txt"
         for source, leg in CROSSING_THRESHOLD_LEG.items() for cat in CROSSING_CATEGORIES],
        # diag_reg_inv_corr: closed-form weighting on the full panel, parallel to
        # free/free_cluster at max-agg only (see compute_diag_reg_inv_corr_weighted_tau)
        expand(f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/figures/tissue_weight_diag_reg_inv_corr.png", source=SOURCES),
        [f"{GENE_SPEC}/{source}/diag_reg_inv_corr_weights/figures/tau_diag_reg_inv_corr_scenarios_scatter_{leg}.png"
         for source, cfg in SOURCES.items() for leg in cfg["legs"]],
        # noise robustness (sigma=0.5): unweighted vs free_cluster_weighted, headline leg/max-agg per source
        expand(f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/figures/tau_noise_variance_compare.png", source=SOURCES),
        # neural genes crossing tau=0.85, GTEx onebrain-max free_cluster
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/crossing_thresholds/onebrain_neural_{NEURAL_CROSSING_THRESHOLD}.csv",
        # tau vs frac_brain 2D histograms, GTEx onebrain-max free_cluster
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/figures/tau_vs_frac_brain_hist2d_onebrain.png",
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/figures/tau_vs_frac_brain_hist2d_onebrain_fracbrain_ge0.4.png",
        # gtex_v11_tx + gtex_v11_tx_protein_coding: transcript-level tau + multi-transcript-gene variance
        expand(f"{GENE_SPEC}/{{txs_source}}/figures/multitranscript_variance.png",
               txs_source=["gtex_v11_tx", "gtex_v11_tx_protein_coding"]),
        expand(f"{GENE_SPEC}/{{txs_source}}/figures/gene_vs_transcript_tau_variance.png",
               txs_source=["gtex_v11_tx", "gtex_v11_tx_protein_coding"]),
        expand(f"{GENE_SPEC}/{{txs_source}}/figures/top_multitranscript_tau_kde.png",
               txs_source=["gtex_v11_tx", "gtex_v11_tx_protein_coding"]),
        expand(f"{GENE_SPEC}/{{txs_source}}/figures/gene_vs_transcript_tau_scatter.png",
               txs_source=["gtex_v11_tx", "gtex_v11_tx_protein_coding"]),
        # gtex_v11_tx: tissue corr mat + tx-vs-gene winning-tissue distance
        f"{GENE_SPEC}/gtex_v11_tx/figures/tx_gene_winner_distance_hist.png",
        # gtex_v11_tx: tissue correlation matrix + its inverse (standard per-source pipeline)
        f"{GENE_SPEC}/gtex_v11_tx/figures/corr_and_inverse.png",
        # gtex_v11_tx: tx-gene winner distance vs. gene/tx tau
        f"{GENE_SPEC}/gtex_v11_tx/figures/winner_distance_vs_tau.png",
        # gtex_v11_tx: gene tau vs. mean transcript tau
        f"{GENE_SPEC}/gtex_v11_tx/figures/gene_tau_vs_mean_tx_tau.png",
        # gtex_v11_tx: log10(TPM+1) vs. tau, gene and tx level, colored by n_transcripts
        f"{GENE_SPEC}/gtex_v11_tx/figures/log_tpm_vs_tau_by_ntx.png",
        # gtex_v11_tx: tx-gene winner distance vs. n_transcripts, colored by gene tau
        f"{GENE_SPEC}/gtex_v11_tx/figures/winner_distance_vs_ntx.png",
        # gtex_v11_tx: DE free_cluster_weights (hierclust0.1, max-agg), tx + gene resolution
        f"{GENE_SPEC}/gtex_v11_tx/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        f"{GENE_SPEC}/gtex_v11_tx_gene/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        # gtex_v11_tx_protein_coding: DE free_cluster_weights (hierclust0.1, max-agg), tx only
        f"{GENE_SPEC}/gtex_v11_tx_protein_coding/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        # gtex_v11_tx_protein_coding: per-gene isoform information content (entropy + KL)
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/effective_n_vs_weighted_kl.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/gene_vs_tx_tau_weighted.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/weighted_kl_vs_tau.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tx_share_hist.png",
        # gtex_v11_tx_protein_coding: donor-level multi-sample bootstrap (SE(tau), SE(TPM), P(tx_tau>gene_tau))
        f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tpm_se.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/gene_tpm_se.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tau_summary.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/gene_tau_summary.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/isoform_fraction_vs_tau.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/gene_tau_by_all_tx_sig.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/specific_isoform_genes.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/specific_isoform_tissue_contribution.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tx_tau_by_gene_tau_decile.png",
        f"{GENE_SPEC}/mxe_chr1_test.csv",
        f"{GENE_SPEC}/transcript_tss.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_vs_tss_variance.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/min_tx_tau_by_gene_tau_decile.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/gene_vs_tx_tau_sig.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tau_deviation_anova_by_tss.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_deviation_anova_by_tss_pval.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_variance_by_tss_count.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_variance_by_tss_count_grouped50bp.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tx_tau_deviation_by_gene_tau_decile.png",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tau_tss_icc.csv",
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tau_tss_icc_arctanh.csv"

# ---- model-fit / distribution-comparison section (unchanged; unrelated to the
# source/gene_filter/agg/method axes above) ----

rule preprocess:
    input:
        GCT
    output:
        f"{PROCESSED}/log10_tpm_long.parquet"
    shell:
        "python3 scripts/preprocess.py {input} {output}"

rule fit_one:
    input:
        f"{PROCESSED}/log10_tpm_long.parquet"
    output:
        f"{MODEL_FITS}/by_model_tissue/{{model}}/{{tissue}}.csv"
    shell:
        "python3 scripts/fit_one.py {wildcards.model} {wildcards.tissue} {input} {output}"

rule combine_model:
    input:
        lambda wc: expand(f"{MODEL_FITS}/by_model_tissue/{wc.model}/{{tissue}}.csv", tissue=TISSUES)
    output:
        f"{MODEL_FITS}/{{model}}_fits.csv"
    shell:
        "python3 scripts/combine_model_fits.py {output} {input}"

rule combine_comparison:
    input:
        expand(f"{MODEL_FITS}/{{model}}_fits.csv", model=MODEL_NAMES)
    output:
        f"{MODEL_FITS}/model_comparison.csv"
    shell:
        "python3 scripts/combine_comparison.py {MODEL_FITS} {output}"

rule plot_model_grid:
    input:
        parquet=f"{PROCESSED}/log10_tpm_long.parquet",
        fits=f"{MODEL_FITS}/{{model}}_fits.csv"
    output:
        f"{DIST_FIGS}/{{model}}_grid.png"
    shell:
        "python3 scripts/plot_model_grid.py {wildcards.model} {input.parquet} {input.fits} {output}"

rule pca_tissue_clustering:
    input:
        f"{MODEL_FITS}/{{model}}_fits.csv"
    output:
        f"{PCA_FIGS}/{{model}}_pca.png"
    shell:
        "python3 scripts/pca_tissue_clustering.py {wildcards.model} {input} {output}"

rule build_expr_matrix:
    input:
        GCT
    output:
        f"{PROCESSED}/log10_tpm_intersection_matrix.parquet"
    shell:
        "python3 scripts/build_expr_matrix.py {input} {output}"

rule expr_vs_model_pca:
    input:
        expr=f"{PROCESSED}/log10_tpm_intersection_matrix.parquet",
        fits=expand(f"{MODEL_FITS}/{{model}}_fits.csv", model=MODEL_NAMES)
    output:
        f"{MODEL_FITS}/pca_compressibility_comparison.csv"
    shell:
        "python3 scripts/expr_vs_model_pca.py {input.expr} {MODEL_FITS} {output}"

rule build_expr_matrix_allgenes:
    input:
        GCT
    output:
        f"{PROCESSED}/log10_tpm_allgenes_floored_matrix.parquet"
    shell:
        "python3 scripts/build_expr_matrix_allgenes.py {input} {output}"

rule expr_vs_model_pca_allgenes:
    input:
        expr=f"{PROCESSED}/log10_tpm_allgenes_floored_matrix.parquet",
        fits=expand(f"{MODEL_FITS}/{{model}}_fits.csv", model=MODEL_NAMES)
    output:
        f"{MODEL_FITS}/pca_compressibility_comparison_allgenes.csv"
    shell:
        "python3 scripts/expr_vs_model_pca.py {input.expr} {MODEL_FITS} {output}"

rule write_methods:
    output:
        "methods.md"
    shell:
        "python3 scripts/write_methods.py {output}"

# ---- reference annotation links + Baron raw-data ingestion (source-specific,
# not wildcarded: genuinely different raw formats) ----

rule link_mane_gtf:
    output:
        f"{DATA}/MANE.v1.5.gtf"
    shell:
        "ln -s {MANE_GTF_SRC} {output}"

rule link_gencode_v47_gtf:
    output:
        GENCODE_V47_GTF_LOCAL
    shell:
        "ln -s {GENCODE_V47_GTF} {output}"

rule build_baron_gene_lengths:
    input:
        mane=f"{DATA}/MANE.v1.5.gtf",
        gencode=GENCODE_V47_GTF,
        baron_csv=BARON_HUMAN_SAMPLES[0]
    output:
        f"{GENE_SPEC_BARON}/gene_lengths.csv"
    shell:
        "python3 scripts/build_baron_gene_lengths.py {input.mane} {input.gencode} {input.baron_csv} {output}"

rule build_baron_celltype_tpm:
    input:
        samples=BARON_HUMAN_SAMPLES,
        lengths=f"{GENE_SPEC_BARON}/gene_lengths.csv"
    output:
        f"{GENE_SPEC_BARON}/celltype_tpm.csv"
    shell:
        "python3 scripts/build_baron_celltype_tpm.py {input.samples} {input.lengths} {output}"

rule build_baron_gct:
    input:
        f"{GENE_SPEC_BARON}/celltype_tpm.csv"
    output:
        f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz"
    shell:
        "python3 scripts/build_baron_gct.py {input} {output}"

# ---- protein-coding filter + full-resolution correlation/clustering (per source) ----

rule filter_protein_coding:
    input:
        gct=lambda wc: RAW_GCT[wc.source],
        gencode=GENCODE_V47_GTF_LOCAL
    output:
        f"{GENE_SPEC}/{{source}}/{{stub}}_protein_coding.gct.gz"
    shell:
        "python3 scripts/filter_protein_coding.py {input.gct} {input.gencode} {output}"

rule compute_tissue_correlation:
    input:
        lambda wc: base_gct(wc.source)
    output:
        corr=f"{GENE_SPEC}/{{source}}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC}/{{source}}/tissue_correlation_linkage.csv"
    shell:
        "python3 scripts/compute_tissue_correlation.py {input} {output.corr} {output.linkage}"

rule plot_tissue_correlation:
    input:
        corr=f"{GENE_SPEC}/{{source}}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC}/{{source}}/tissue_correlation_linkage.csv"
    output:
        f"{GENE_SPEC}/{{source}}/figures/tissue_correlation_clustermap.png"
    params:
        title=lambda wc: TAU_PLOT_TITLE[wc.source]
    shell:
        """python3 scripts/plot_tissue_correlation.py {input.corr} {input.linkage} {output} "{params.title}" """

rule plot_corr_and_inverse:
    input:
        corr=f"{GENE_SPEC}/{{source}}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC}/{{source}}/tissue_correlation_linkage.csv"
    output:
        f"{GENE_SPEC}/{{source}}/figures/corr_and_inverse.png"
    params:
        title=lambda wc: TAU_PLOT_TITLE[wc.source]
    shell:
        """python3 scripts/plot_corr_and_inverse.py {input.corr} {input.linkage} {output} "{params.title}" """

rule compute_tau_cluster_weighted_full_gtex:
    input:
        gct=lambda wc: base_gct("gtex_protein_coding"),
        linkage=f"{GENE_SPEC_GTEX_PC}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_GTEX_PC}/tissue_correlation.csv"
    output:
        clusters=CLUSTERS_CSV["gtex_protein_coding"],
        tau=CLUSTERS_TAU_CSV["gtex_protein_coding"]
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} {output.clusters} {output.tau}"

rule compute_tau_cluster_weighted_full_baron:
    input:
        gct=lambda wc: base_gct("baron_protein_coding"),
        linkage=f"{GENE_SPEC_BARON_PC}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_BARON_PC}/tissue_correlation.csv"
    output:
        clusters=CLUSTERS_CSV["baron_protein_coding"],
        tau=CLUSTERS_TAU_CSV["baron_protein_coding"]
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} {output.clusters} {output.tau}"

rule plot_tissue_clusters_dendrogram:
    input:
        linkage=f"{GENE_SPEC}/{{source}}/tissue_correlation_linkage.csv",
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        f"{GENE_SPEC}/{{source}}/figures/tissue_clusters_dendrogram_hierclust0.1.png"
    params:
        label=lambda wc: SOURCES[wc.source]["display"]
    shell:
        "python3 scripts/plot_tissue_clusters_dendrogram.py {input.linkage} {input.clusters} {output} '{params.label}'"

# ---- leg (resolution) GCTs + reference weighted_tau, wildcarded over source/leg/agg ----

rule build_leg_gct:
    input:
        gct=lambda wc: base_gct(wc.source)
    output:
        f"{GENE_SPEC}/{{source}}/{{leg}}{{agg_suffix}}/{{stub}}_{{leg}}.gct.gz"
    wildcard_constraints:
        leg=r"onebrain|oneislets"
    params:
        agg=lambda wc: wc.agg_suffix.lstrip("_"),
        members=lambda wc: " ".join(leg_members(wc.source, wc.leg))
    shell:
        "python3 scripts/build_gct_onegroup.py {input.gct} {output} {wildcards.leg} {params.agg} {params.members}"

rule build_cluster_leg_gct:
    input:
        gct=lambda wc: base_gct(wc.source),
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        f"{GENE_SPEC}/{{source}}/hierclust0.1{{agg_suffix}}/{{stub}}_hierclust0.1.gct.gz"
    params:
        agg=lambda wc: wc.agg_suffix.lstrip("_")
    shell:
        "python3 scripts/build_gct_clustered.py {input.gct} {input.clusters} {output} {params.agg}"

rule compute_weighted_tau_leg:
    input:
        lambda wc: leg_gct_path(wc.source, wc.leg, wc.agg_suffix)
    output:
        tau=f"{GENE_SPEC}/{{source}}/{{leg}}{{agg_suffix}}/weighted_tau.csv",
        weights=f"{GENE_SPEC}/{{source}}/{{leg}}{{agg_suffix}}/weighted_tau_weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau.py {input} {output.tau} {output.weights}"

# ---- free_weights / free_cluster_weights: DE chains, wildcarded over
# source/leg/agg/chain. Two methods kept as separate rules (not a further
# {method} wildcard) since free_cluster_weights has one extra input (clusters)
# and one extra output (cluster_weights.csv) that free_weights doesn't. ----

rule optimize_free_weights_chain:
    input:
        full_gct=lambda wc: base_gct(wc.source),
        ref_tau=lambda wc: leg_tau_path(wc.source, wc.leg, wc.agg_suffix)
    output:
        weights=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/weights.csv",
        tau=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/weighted_tau.csv",
        log=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/optimization_log.csv",
        best=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/best.csv"
    shell:
        "python3 scripts/optimize_free_weights.py tissue {input.full_gct} {input.ref_tau} - "
        "{output.weights} - {output.tau} {output.log} {output.best} {wildcards.chain}"

rule optimize_free_cluster_weights_chain:
    input:
        full_gct=lambda wc: base_gct(wc.source),
        ref_tau=lambda wc: leg_tau_path(wc.source, wc.leg, wc.agg_suffix),
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        weights=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/weights.csv",
        cluster_weights=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/cluster_weights.csv",
        tau=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/weighted_tau.csv",
        log=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/optimization_log.csv",
        best=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/chains/chain_{{chain}}/best.csv"
    shell:
        "python3 scripts/optimize_free_weights.py cluster {input.full_gct} {input.ref_tau} {input.clusters} "
        "{output.weights} {output.cluster_weights} {output.tau} {output.log} {output.best} {wildcards.chain}"

rule combine_free_weights_chains:
    input:
        lambda wc: expand(f"{GENE_SPEC}/{wc.source}/free_weights{wc.agg_suffix}/{wc.leg}/chains/chain_{{chain}}/best.csv",
                           chain=range(SOURCES[wc.source]["n_free"]))
    output:
        weights=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/weights.csv",
        tau=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/weighted_tau.csv",
        summary=f"{GENE_SPEC}/{{source}}/free_weights{{agg_suffix}}/{{leg}}/chains_summary.csv"
    params:
        chains_dir=lambda wc: f"{GENE_SPEC}/{wc.source}/free_weights{wc.agg_suffix}/{wc.leg}/chains",
        n=lambda wc: SOURCES[wc.source]["n_free"]
    shell:
        "python3 scripts/combine_free_weights_chains.py tissue {params.chains_dir} {params.n} "
        "{output.weights} - {output.tau} {output.summary}"

rule combine_free_cluster_weights_chains:
    input:
        lambda wc: expand(f"{GENE_SPEC}/{wc.source}/free_cluster_weights{wc.agg_suffix}/{wc.leg}/chains/chain_{{chain}}/best.csv",
                           chain=range(SOURCES[wc.source]["n_free_cluster"]))
    output:
        weights=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/weights.csv",
        cluster_weights=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/cluster_weights.csv",
        tau=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/weighted_tau.csv",
        summary=f"{GENE_SPEC}/{{source}}/free_cluster_weights{{agg_suffix}}/{{leg}}/chains_summary.csv"
    params:
        chains_dir=lambda wc: f"{GENE_SPEC}/{wc.source}/free_cluster_weights{wc.agg_suffix}/{wc.leg}/chains",
        n=lambda wc: SOURCES[wc.source]["n_free_cluster"]
    shell:
        "python3 scripts/combine_free_weights_chains.py cluster {params.chains_dir} {params.n} "
        "{output.weights} {output.cluster_weights} {output.tau} {output.summary}"

# ---- plotting: weights bar chart (source-generic script, one rule -- mean and
# max nest identically under {method}{agg_suffix}/figures/) ----

rule plot_free_weights_clean:
    input:
        weights=lambda wc: f"{GENE_SPEC}/{wc.source}/free{wc.cluster_tag}_weights{wc.agg_suffix}/{wc.leg}/weights.csv",
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        f"{GENE_SPEC}/{{source}}/free{{cluster_tag}}_weights{{agg_suffix}}/figures/tissue_weight_free{{cluster_tag}}_{{leg}}.png"
    wildcard_constraints:
        cluster_tag=r"|_cluster"
    shell:
        "python3 scripts/plot_free_weights.py {input.weights} {input.clusters} {output}"

# ---- plotting: tau-vs-reference comparison. GTEx and Baron use different
# scripts (plot_tau_free_vs_full.py needs the full-resolution GCT to also show
# tau_unweighted on both sides; plot_tau_scenarios_scatter.py doesn't) -- kept
# as separate small rules per source, each still wildcarded over leg/agg. ----

rule plot_tau_free_vs_full_gtex:
    input:
        gct=lambda wc: base_gct("gtex_protein_coding"),
        free_tau=f"{GENE_SPEC_GTEX_PC}/free_weights{{agg_suffix}}/{{leg}}/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path("gtex_protein_coding", wc.leg, wc.agg_suffix)
    output:
        f"{GENE_SPEC_GTEX_PC}/free_weights{{agg_suffix}}/figures/tau_free_vs_full_{{leg}}.png"
    shell:
        "python3 scripts/plot_tau_free_vs_full.py {input.gct} {input.free_tau} {input.ref_tau} {output} {wildcards.leg}"

rule plot_tau_free_cluster_vs_full_gtex:
    input:
        gct=lambda wc: base_gct("gtex_protein_coding"),
        free_tau=f"{GENE_SPEC_GTEX_PC}/free_cluster_weights{{agg_suffix}}/{{leg}}/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path("gtex_protein_coding", wc.leg, wc.agg_suffix)
    output:
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights{{agg_suffix}}/figures/tau_free_cluster_vs_full_{{leg}}.png"
    shell:
        "python3 scripts/plot_tau_free_vs_full.py {input.gct} {input.free_tau} {input.ref_tau} {output} {wildcards.leg}"

rule plot_tau_scenarios_baron:
    input:
        free_tau=f"{GENE_SPEC_BARON_PC}/free_weights{{agg_suffix}}/{{leg}}/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path("baron_protein_coding", wc.leg, wc.agg_suffix)
    output:
        f"{GENE_SPEC_BARON_PC}/free_weights{{agg_suffix}}/figures/tau_scenarios_scatter_{{leg}}.png"
    shell:
        "python3 scripts/plot_tau_scenarios_scatter.py {input.free_tau} {input.ref_tau} {output} {wildcards.leg}"

rule plot_tau_free_cluster_scenarios_baron:
    input:
        free_tau=f"{GENE_SPEC_BARON_PC}/free_cluster_weights{{agg_suffix}}/{{leg}}/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path("baron_protein_coding", wc.leg, wc.agg_suffix)
    output:
        f"{GENE_SPEC_BARON_PC}/free_cluster_weights{{agg_suffix}}/figures/tau_free_cluster_scenarios_scatter_{{leg}}.png"
    shell:
        "python3 scripts/plot_tau_scenarios_scatter.py {input.free_tau} {input.ref_tau} {output} {wildcards.leg}"

rule list_crossing_threshold_genes:
    input:
        f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/{{leg}}/weighted_tau.csv"
    output:
        [f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/crossing_thresholds/{{leg}}_{cat}.txt"
         for cat in CROSSING_CATEGORIES]
    params:
        out_dir=f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/crossing_thresholds"
    wildcard_constraints:
        source=r"gtex_protein_coding|baron_protein_coding",
        leg=r"onebrain|oneislets"
    shell:
        "python3 scripts/list_crossing_threshold_genes.py {input} {params.out_dir} {wildcards.leg}"

# ---- diag_reg_inv_corr ("diagonal-regularized inverse correlation"): closed-form
# weighting (compute_weighted_tau.py, already used elsewhere for leg-local reference
# tau) applied to the full base_gct, parallel to free_weights/free_cluster_weights but
# agg/leg-independent -- one weight vector per source, reused across every leg's
# reference in the scenario plot. Max-agg only for now (see
# plot_tau_diag_reg_inv_corr_scenarios). ----

rule compute_diag_reg_inv_corr_weighted_tau:
    input:
        lambda wc: base_gct(wc.source)
    output:
        tau=f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/weighted_tau.csv",
        weights=f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau.py {input} {output.tau} {output.weights}"

rule plot_diag_reg_inv_corr_weights:
    input:
        weights=f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/weights.csv",
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/figures/tissue_weight_diag_reg_inv_corr.png"
    shell:
        "python3 scripts/plot_diag_reg_inv_corr_weights.py {input.weights} {input.clusters} {output}"

rule plot_tau_diag_reg_inv_corr_scenarios:
    input:
        diag_reg_tau=f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path(wc.source, wc.leg, "_max")
    output:
        f"{GENE_SPEC}/{{source}}/diag_reg_inv_corr_weights/figures/tau_diag_reg_inv_corr_scenarios_scatter_{{leg}}.png"
    shell:
        "python3 scripts/plot_tau_scenarios_scatter.py {input.diag_reg_tau} {input.ref_tau} {output} {wildcards.leg} "
        "tau_weighted_diag_reg_inv_corr diagonal-regularized-inverse-corr"

# ---- margin: closed-form weighting directly from the expression matrix (no
# correlation matrix at all, see gls_tau.margin_weights) applied to the full
# base_gct, parallel to diag_reg_inv_corr/free/free_cluster but agg/leg-independent
# and reference-free. Max-agg only for now, same as diag_reg_inv_corr. ----

rule compute_margin_weighted_tau:
    input:
        lambda wc: base_gct(wc.source)
    output:
        tau=f"{GENE_SPEC}/{{source}}/margin_weights/weighted_tau.csv",
        weights=f"{GENE_SPEC}/{{source}}/margin_weights/weights.csv"
    shell:
        "python3 scripts/compute_margin_weighted_tau.py {input} {output.tau} {output.weights}"

rule plot_margin_weights:
    input:
        weights=f"{GENE_SPEC}/{{source}}/margin_weights/weights.csv",
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        f"{GENE_SPEC}/{{source}}/margin_weights/figures/tissue_weight_margin.png"
    shell:
        "python3 scripts/plot_margin_weights.py {input.weights} {input.clusters} {output}"

rule plot_tau_margin_scenarios:
    input:
        margin_tau=f"{GENE_SPEC}/{{source}}/margin_weights/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path(wc.source, wc.leg, "_max")
    output:
        f"{GENE_SPEC}/{{source}}/margin_weights/figures/tau_margin_scenarios_scatter_{{leg}}.png"
    shell:
        "python3 scripts/plot_tau_scenarios_scatter.py {input.margin_tau} {input.ref_tau} {output} {wildcards.leg} "
        "tau_weighted_margin margin-weighted"

# ---- noise robustness: Monte Carlo perturbation of log10(TPM+1) (see
# gls_tau.apply_nonzero_log_noise, which clips noisy TPM at 0), weights frozen at
# their already-fitted value (not refit per replicate -- see applying_noise.md).
# free_cluster_weights_max at each source's headline leg (CROSSING_THRESHOLD_LEG),
# sigma=0.5 (an earlier sigma=0.1 pass was retired in favor of this higher-sigma
# run). ----

rule compute_noise_variance_free_cluster:
    input:
        gct=lambda wc: base_gct(wc.source),
        weights=lambda wc: f"{GENE_SPEC}/{wc.source}/free_cluster_weights_max/{CROSSING_THRESHOLD_LEG[wc.source]}/weights.csv"
    output:
        summary=f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/noise_variance/summary.csv",
        long=f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/noise_variance/pergene_long.csv"
    params:
        n=NOISE_N_REPLICATES,
        sigma=NOISE_SIGMA,
        seed=NOISE_SEED
    shell:
        "python3 scripts/compute_noise_variance.py {input.gct} {input.weights} weight free_cluster "
        "{params.n} {params.sigma} {params.seed} {output.summary} {output.long}"

rule plot_noise_variance_compare_free_cluster:
    input:
        f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/noise_variance/summary.csv"
    output:
        f"{GENE_SPEC}/{{source}}/free_cluster_weights_max/figures/tau_noise_variance_compare.png"
    params:
        label=lambda wc: SOURCES[wc.source]["display"]
    shell:
        "python3 scripts/plot_jackknife_variance_compare.py {input} {output} '{params.label}' "
        "tau_unweighted_noise_var tau_free_cluster_weighted_noise_var noise free_cluster_weighted"

# ---- GENCODE gene_id -> gene_name reference map + neural genes crossing tau=0.85
# under free_cluster weighting (GTEx onebrain-max headline leg). ----

rule build_gencode_gene_id_name_map:
    input:
        GENCODE_V47_GTF_LOCAL
    output:
        f"{GENE_SPEC}/gencode_v47_gene_id_name_map.csv"
    shell:
        "python3 scripts/build_gencode_gene_id_name_map.py {input} {output}"

rule extract_transcript_tss:
    input:
        GENCODE_V47_GTF_LOCAL
    output:
        f"{GENE_SPEC}/transcript_tss.csv"
    shell:
        "python3 scripts/extract_transcript_tss.py {input} {output}"

rule list_neural_threshold_crossing_genes:
    input:
        gct=lambda wc: base_gct("gtex_protein_coding"),
        tau=f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/onebrain/weighted_tau.csv",
        name_map=f"{GENE_SPEC}/gencode_v47_gene_id_name_map.csv"
    output:
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/crossing_thresholds/onebrain_neural_{NEURAL_CROSSING_THRESHOLD}.csv"
    shell:
        "python3 scripts/list_neural_threshold_crossing_genes.py {input.gct} {input.tau} {input.name_map} "
        f"{NEURAL_CROSSING_THRESHOLD} " "{output}"

rule plot_tau_vs_frac_brain_hist2d:
    input:
        gct=lambda wc: base_gct("gtex_protein_coding"),
        tau=f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/onebrain/weighted_tau.csv"
    output:
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/figures/tau_vs_frac_brain_hist2d_onebrain.png"
    shell:
        "python3 scripts/plot_tau_vs_frac_brain_hist2d.py {input.gct} {input.tau} {output}"

rule plot_tau_vs_frac_brain_hist2d_minfrac:
    input:
        gct=lambda wc: base_gct("gtex_protein_coding"),
        tau=f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/onebrain/weighted_tau.csv"
    output:
        f"{GENE_SPEC_GTEX_PC}/free_cluster_weights_max/figures/tau_vs_frac_brain_hist2d_onebrain_fracbrain_ge{{min_frac}}.png"
    shell:
        "python3 scripts/plot_tau_vs_frac_brain_hist2d.py {input.gct} {input.tau} {output} {wildcards.min_frac}"

# ---- gtex_v11_tx raw-data ingestion (RSEM quant, full 68-tissue panel, Ensembl IDs) ----

rule download_gtex_v11_tx_tpm:
    output:
        GTEX_V11_TX_TPM_RAW
    shell:
        "wget -O {output} {GTEX_V11_TX_TPM_URL}"

rule download_gtex_v11_sample_attrs:
    output:
        GTEX_V11_SAMPLE_ATTRS_RAW
    shell:
        "wget -O {output} {GTEX_V11_SAMPLE_ATTRS_URL}"

rule build_gtex_v11_sample_tissue_map:
    input:
        GTEX_V11_SAMPLE_ATTRS_RAW
    output:
        f"{GENE_SPEC_GTEX_V11_TX}/sample_tissue_map.csv"
    shell:
        "python3 scripts/build_gtex_v11_sample_tissue_map.py {input} {output}"

rule build_gtex_v11_tx_gct:
    input:
        tpm=GTEX_V11_TX_TPM_RAW,
        sample_tissue=f"{GENE_SPEC_GTEX_V11_TX}/sample_tissue_map.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX}/transcript_median_tpm.gct.gz"
    shell:
        "python3 scripts/build_gtex_v11_tx_gct.py {input.tpm} {input.sample_tissue} {output}"

rule filter_protein_coding_transcripts_gtex_v11_tx:
    input:
        gct=f"{GENE_SPEC_GTEX_V11_TX}/transcript_median_tpm.gct.gz",
        gencode=GENCODE_V47_GTF_LOCAL
    output:
        GTEX_V11_TX_BASE_GCT
    shell:
        "python3 scripts/filter_protein_coding_by_gene_id.py {input.gct} {input.gencode} {output}"

rule build_gtex_v11_tx_meta:
    input:
        gct=GTEX_V11_TX_BASE_GCT,
        name_map=f"{GENE_SPEC}/gencode_v47_gene_id_name_map.csv"
    output:
        TXS_META_CSV["gtex_v11_tx"]
    shell:
        "python3 scripts/build_gtex_v11_tx_meta.py {input.gct} {input.name_map} {output}"

rule filter_protein_coding_transcripts_gtex_v11_tx_pc:
    input:
        gct=GTEX_V11_TX_BASE_GCT,
        gencode=GENCODE_V47_GTF_LOCAL
    output:
        GTEX_V11_TX_PC_BASE_GCT
    shell:
        "python3 scripts/filter_protein_coding_transcript_biotype.py {input.gct} {input.gencode} {output}"

rule build_gtex_v11_tx_pc_meta:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        name_map=f"{GENE_SPEC}/gencode_v47_gene_id_name_map.csv"
    output:
        TXS_META_CSV["gtex_v11_tx_protein_coding"]
    shell:
        "python3 scripts/build_gtex_v11_tx_meta.py {input.gct} {input.name_map} {output}"

rule compute_tx_info_content:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        weights=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weights.csv",
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv"
    output:
        pickle=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_kl_divergence.pkl",
        csv=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    shell:
        "python3 scripts/compute_tx_info_content.py {input.gct} {input.weights} {input.tx_tau} {output.pickle} {output.csv}"

rule plot_effective_n_vs_weighted_kl:
    input:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/effective_n_vs_weighted_kl.png"
    shell:
        "python3 scripts/plot_effective_n_vs_weighted_kl.py {input} {output}"

rule plot_gene_vs_tx_tau_weighted:
    input:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/gene_vs_tx_tau_weighted.png"
    shell:
        "python3 scripts/plot_gene_vs_tx_tau_weighted.py {input} {output}"

rule plot_weighted_kl_vs_tau:
    input:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/weighted_kl_vs_tau.png"
    shell:
        "python3 scripts/plot_weighted_kl_vs_tau.py {input} {output}"

rule plot_tx_share_hist:
    input:
        GTEX_V11_TX_PC_BASE_GCT
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tx_share_hist.png"
    shell:
        "python3 scripts/plot_tx_share_hist.py {input} {output}"

N_BOOTSTRAP_GTEX_V11_TX_PC = 200
BOOTSTRAP_SEED_GTEX_V11_TX_PC = 0

rule plot_isoform_fraction_vs_tau:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv",
        bootstrap_tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tau_summary.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/isoform_fraction_vs_tau.png"
    shell:
        "python3 scripts/plot_isoform_fraction_vs_tau.py {input.gct} {input.tx_tau} {input.info_content} "
        "{input.bootstrap_tx_tau} {output}"

rule plot_gene_tau_by_all_tx_sig:
    input:
        bootstrap_tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tau_summary.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/gene_tau_by_all_tx_sig.png"
    shell:
        "python3 scripts/plot_gene_tau_by_all_tx_sig.py {input.bootstrap_tx_tau} {input.info_content} {output}"

rule find_specific_isoform_genes:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        bootstrap_tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tau_summary.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/specific_isoform_genes.csv"
    shell:
        "python3 scripts/find_specific_isoform_genes.py {input.gct} {input.bootstrap_tx_tau} {input.info_content} {output}"

rule plot_specific_isoform_tissue_contribution:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        specific_isoform=f"{GENE_SPEC_GTEX_V11_TX_PC}/specific_isoform_genes.csv",
        clusters=CLUSTERS_CSV["gtex_v11_tx_protein_coding"]
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/specific_isoform_tissue_contribution.png"
    shell:
        "python3 scripts/plot_specific_isoform_tissue_contribution.py {input.gct} {input.specific_isoform} "
        "{input.clusters} {output}"

rule plot_tx_tau_by_gene_tau_decile:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tx_tau_by_gene_tau_decile.png"
    shell:
        "python3 scripts/plot_tx_tau_by_gene_tau_decile.py {input.gct} {input.tx_tau} {input.info_content} {output}"

rule plot_min_tx_tau_by_gene_tau_decile:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/min_tx_tau_by_gene_tau_decile.png"
    shell:
        "python3 scripts/plot_min_tx_tau_by_gene_tau_decile.py {input.gct} {input.tx_tau} {input.info_content} {output}"

rule plot_gene_vs_tx_tau_sig:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/gene_vs_tx_tau_sig.png"
    shell:
        "python3 scripts/plot_gene_vs_tx_tau_sig.py {input.gct} {input.tx_tau} {input.info_content} {output}"

rule anova_tau_deviation_by_tss:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv",
        tss=f"{GENE_SPEC}/transcript_tss.csv"
    output:
        csv=f"{GENE_SPEC_GTEX_V11_TX_PC}/tau_deviation_anova_by_tss.csv",
        png=f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_deviation_anova_by_tss_pval.png"
    shell:
        "python3 scripts/anova_tau_deviation_by_tss.py {input.gct} {input.tx_tau} {input.info_content} {input.tss} "
        "{output.csv} {output.png}"

rule plot_tau_variance_by_tss_count:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        tss=f"{GENE_SPEC}/transcript_tss.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_variance_by_tss_count.png"
    shell:
        "python3 scripts/plot_tau_variance_by_tss_count.py {input.gct} {input.tx_tau} {input.tss} {output}"

rule plot_tau_variance_by_tss_count_grouped50bp:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        tss=f"{GENE_SPEC}/transcript_tss.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_variance_by_tss_count_grouped50bp.png"
    shell:
        "python3 scripts/plot_tau_variance_by_tss_count_grouped50bp.py {input.gct} {input.tx_tau} {input.tss} {output}"

rule plot_tx_tau_deviation_by_gene_tau_decile:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tx_tau_deviation_by_gene_tau_decile.png"
    shell:
        "python3 scripts/plot_tx_tau_deviation_by_gene_tau_decile.py {input.gct} {input.tx_tau} {input.info_content} {output}"

rule compute_tau_tss_icc:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv",
        tss=f"{GENE_SPEC}/transcript_tss.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tau_tss_icc.csv"
    shell:
        "python3 scripts/compute_tau_tss_icc.py {input.gct} {input.tx_tau} {input.info_content} {input.tss} {output}"

rule compute_tau_tss_icc_arctanh:
    input:
        gct=GTEX_V11_TX_PC_BASE_GCT,
        tx_tau=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv",
        tss=f"{GENE_SPEC}/transcript_tss.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/tau_tss_icc_arctanh.csv"
    shell:
        "python3 scripts/compute_tau_tss_icc_arctanh.py {input.gct} {input.tx_tau} {input.info_content} {input.tss} {output}"

rule find_mutually_exclusive_exons_chr1:
    input:
        GENCODE_V47_GTF_LOCAL
    output:
        f"{GENE_SPEC}/mxe_chr1_test.csv"
    shell:
        "python3 scripts/find_mutually_exclusive_exons.py {input} {output} chr1"

rule plot_tau_vs_tss_variance:
    input:
        tss=f"{GENE_SPEC}/transcript_tss.csv",
        transcript_meta=TXS_META_CSV["gtex_v11_tx_protein_coding"],
        info_content=f"{GENE_SPEC_GTEX_V11_TX_PC}/tx_info_content.csv",
        specific_isoform=f"{GENE_SPEC_GTEX_V11_TX_PC}/specific_isoform_genes.csv"
    output:
        f"{GENE_SPEC_GTEX_V11_TX_PC}/figures/tau_vs_tss_variance.png"
    shell:
        "python3 scripts/plot_tau_vs_tss_variance.py {input.tss} {input.transcript_meta} {input.info_content} "
        "{input.specific_isoform} {output}"

rule bootstrap_gtex_v11_tx_protein_coding:
    input:
        raw_tpm=GTEX_V11_TX_TPM_RAW,
        sample_tissue=f"{GENE_SPEC_GTEX_V11_TX}/sample_tissue_map.csv",
        gct=GTEX_V11_TX_PC_BASE_GCT,
        weights=f"{GENE_SPEC_GTEX_V11_TX_PC}/free_cluster_weights_max/hierclust0.1/weights.csv"
    output:
        tx_tpm_se=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tpm_se.csv",
        gene_tpm_se=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/gene_tpm_se.csv",
        tx_tau_summary=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/tx_tau_summary.csv",
        gene_tau_summary=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap/gene_tau_summary.csv"
    params:
        out_dir=f"{GENE_SPEC_GTEX_V11_TX_PC}/bootstrap",
        n=N_BOOTSTRAP_GTEX_V11_TX_PC,
        seed=BOOTSTRAP_SEED_GTEX_V11_TX_PC
    shell:
        "python3 scripts/bootstrap_gtex_v11_tx_protein_coding.py {input.raw_tpm} {input.sample_tissue} "
        "{input.gct} {input.weights} {params.out_dir} {params.n} {params.seed}"

# ---- shared downstream analysis (tau / multi-transcript variance / figures),
# identical in shape for both sources -- wildcarded over txs_source. ----

rule compute_tau_unweighted_txs:
    input:
        lambda wc: TXS_BASE_GCT[wc.txs_source]
    output:
        f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/compute_tau_unweighted.py {input} {output}"

rule compute_multitranscript_variance:
    input:
        gct=lambda wc: TXS_BASE_GCT[wc.txs_source],
        tau=f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/multitranscript_variance.csv"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/compute_multitranscript_variance.py {input.gct} {input.tau} {output}"

rule plot_multitranscript_variance:
    input:
        f"{GENE_SPEC}/{{txs_source}}/multitranscript_variance.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/multitranscript_variance.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_multitranscript_variance.py {input} {output}"

rule compute_gene_tau_txs:
    input:
        lambda wc: TXS_BASE_GCT[wc.txs_source]
    output:
        f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/compute_gene_tau_gtex_txs.py {input} {output}"

rule plot_gene_vs_transcript_tau_variance:
    input:
        gene_tau=f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv",
        multitx=f"{GENE_SPEC}/{{txs_source}}/multitranscript_variance.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/gene_vs_transcript_tau_variance.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_gene_vs_transcript_tau_variance.py {input.gene_tau} {input.multitx} {output}"

rule plot_top_multitranscript_tau_kde:
    input:
        tau=f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv",
        meta=lambda wc: TXS_META_CSV[wc.txs_source],
        multitx=f"{GENE_SPEC}/{{txs_source}}/multitranscript_variance.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/top_multitranscript_tau_kde.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_top_multitranscript_tau_kde.py {input.tau} {input.meta} {input.multitx} {output}"

rule plot_gene_vs_transcript_tau_scatter:
    input:
        tx_tau=f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv",
        meta=lambda wc: TXS_META_CSV[wc.txs_source],
        gene_tau=f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv",
        multitx=f"{GENE_SPEC}/{{txs_source}}/multitranscript_variance.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/gene_vs_transcript_tau_scatter.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_gene_vs_transcript_tau_scatter.py {input.tx_tau} {input.meta} {input.gene_tau} {input.multitx} {output} {wildcards.txs_source}"

rule compute_tx_gene_winner_distance:
    input:
        lambda wc: TXS_BASE_GCT[wc.txs_source]
    output:
        corr=f"{GENE_SPEC}/{{txs_source}}/tissue_correlation_matrix.csv",
        dist=f"{GENE_SPEC}/{{txs_source}}/tx_gene_winner_distance.csv"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/compute_tx_gene_winner_distance.py {input} {output.corr} {output.dist}"

rule plot_tx_gene_winner_distance_hist:
    input:
        f"{GENE_SPEC}/{{txs_source}}/tx_gene_winner_distance.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/tx_gene_winner_distance_hist.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_tx_gene_winner_distance_hist.py {input} {output}"

rule plot_winner_distance_vs_tau:
    input:
        dist=f"{GENE_SPEC}/{{txs_source}}/tx_gene_winner_distance.csv",
        gene_tau=f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv",
        tx_tau=f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv"
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/winner_distance_vs_tau.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_winner_distance_vs_tau.py {input.dist} {input.gene_tau} {input.tx_tau} {output} {wildcards.txs_source}"

rule plot_gene_tau_vs_mean_tx_tau:
    input:
        tx_tau=f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv",
        meta=lambda wc: TXS_META_CSV[wc.txs_source],
        gene_tau=f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv",
        multitx=f"{GENE_SPEC}/{{txs_source}}/multitranscript_variance.csv",
        gct=lambda wc: TXS_BASE_GCT[wc.txs_source]
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/gene_tau_vs_mean_tx_tau.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_gene_tau_vs_mean_tx_tau.py {input.tx_tau} {input.meta} {input.gene_tau} {input.multitx} {input.gct} {output} {wildcards.txs_source}"

rule plot_log_tpm_vs_tau_by_ntx:
    input:
        gene_tau=f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv",
        tx_tau=f"{GENE_SPEC}/{{txs_source}}/weighted_tau.csv",
        gct=lambda wc: TXS_BASE_GCT[wc.txs_source]
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/log_tpm_vs_tau_by_ntx.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_log_tpm_vs_tau_by_ntx.py {input.gene_tau} {input.tx_tau} {input.gct} {output} {wildcards.txs_source}"

rule plot_winner_distance_vs_ntx:
    input:
        dist=f"{GENE_SPEC}/{{txs_source}}/tx_gene_winner_distance.csv",
        gene_tau=f"{GENE_SPEC}/{{txs_source}}/gene_tau.csv",
        gct=lambda wc: TXS_BASE_GCT[wc.txs_source]
    output:
        f"{GENE_SPEC}/{{txs_source}}/figures/winner_distance_vs_ntx.png"
    wildcard_constraints:
        txs_source = TXS_SOURCES_RE
    shell:
        "python3 scripts/plot_winner_distance_vs_ntx.py {input.dist} {input.gene_tau} {input.gct} {output} {wildcards.txs_source}"

# ---- gtex_v11_tx DE free_cluster_weights (tx + gene resolution) -- see block
# comment near GTEX_V11_TX_GENE_BASE_GCT for design rationale ----

rule build_gtex_v11_tx_gene_gct:
    input:
        GTEX_V11_TX_BASE_GCT
    output:
        GTEX_V11_TX_GENE_BASE_GCT
    shell:
        "python3 scripts/build_gtex_v11_tx_gene_gct.py {input} {output}"

rule compute_tau_cluster_weighted_gtex_v11_tx:
    input:
        gct=lambda wc: base_gct(wc.res),
        linkage=lambda wc: f"{GENE_SPEC}/{CLUSTER_CORR_BASIS[wc.res]}/tissue_correlation_linkage.csv",
        corr=lambda wc: f"{GENE_SPEC}/{CLUSTER_CORR_BASIS[wc.res]}/tissue_correlation.csv"
    output:
        clusters=f"{GENE_SPEC}/{{res}}/cluster_weights/tissue_clusters_hierclust0.1.csv",
        tau=f"{GENE_SPEC}/{{res}}/cluster_weights/weighted_tau_hierclust0.1.csv"
    wildcard_constraints:
        res = GTEX_V11_TX_RES_RE
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} {output.clusters} {output.tau}"

rule combine_free_cluster_weights_chains_gtex_v11_tx:
    input:
        lambda wc: expand(f"{GENE_SPEC}/{wc.res}/free_cluster_weights_max/hierclust0.1/chains/chain_{{chain}}/best.csv",
                           chain=range(N_FREE_CLUSTER_GTEX_V11_TX))
    output:
        weights=f"{GENE_SPEC}/{{res}}/free_cluster_weights_max/hierclust0.1/weights.csv",
        cluster_weights=f"{GENE_SPEC}/{{res}}/free_cluster_weights_max/hierclust0.1/cluster_weights.csv",
        tau=f"{GENE_SPEC}/{{res}}/free_cluster_weights_max/hierclust0.1/weighted_tau.csv",
        summary=f"{GENE_SPEC}/{{res}}/free_cluster_weights_max/hierclust0.1/chains_summary.csv"
    wildcard_constraints:
        res = GTEX_V11_TX_RES_RE
    params:
        chains_dir=lambda wc: f"{GENE_SPEC}/{wc.res}/free_cluster_weights_max/hierclust0.1/chains"
    shell:
        "python3 scripts/combine_free_weights_chains.py cluster {params.chains_dir} {N_FREE_CLUSTER_GTEX_V11_TX} "
        "{output.weights} {output.cluster_weights} {output.tau} {output.summary}"
