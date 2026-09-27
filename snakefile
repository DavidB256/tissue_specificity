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
#   resolution/"leg" (onebrain|oneislets | simplified_0.1) -- per-source, see
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
        legs={"onebrain": ("prefix", "Brain"), "simplified_0.1": ("cluster", None)},
        n_free=60, n_free_cluster=20, display="GTEx (protein-coding)",
    ),
    "baron_protein_coding": dict(
        legs={"oneislets": ("group", BARON_ONEISLETS_MEMBERS), "simplified_0.1": ("cluster", None)},
        n_free=10, n_free_cluster=10, display="Baron pancreas (protein-coding)",
    ),
}
STUB = {"gtex_protein_coding": "gene_median_tpm", "baron_protein_coding": "celltype_tpm"}
RAW_GCT = {"gtex_protein_coding": GCT, "baron_protein_coding": f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz"}
# legacy paths (predate this wildcard scheme, kept as literals rather than migrated --
# no data movement, see streamline-snakefile.md plan): GTEx nests clusters under a
# "full/" subdir, Baron doesn't.
CLUSTERS_CSV = {
    "gtex_protein_coding": f"{GENE_SPEC_GTEX_PC}/cluster_weights/full/tissue_clusters.csv",
    "baron_protein_coding": f"{GENE_SPEC_BARON_PC}/cluster_weights/tissue_clusters_simplified_0.1.csv",
}
CLUSTERS_TAU_CSV = {
    "gtex_protein_coding": f"{GENE_SPEC_GTEX_PC}/cluster_weights/full/weighted_tau.csv",
    "baron_protein_coding": f"{GENE_SPEC_BARON_PC}/cluster_weights/weighted_tau_simplified_0.1.csv",
}
TAU_PLOT_TITLE = {
    "gtex_protein_coding": "log10 TPM, TPM<1->0, GTEx tissues (protein-coding)",
    "baron_protein_coding": "log10 TPM, TPM<1->0, Baron pancreas cell types (protein-coding)",
}
# headline free_cluster_weights_max leg per source, for the crossing_thresholds gene lists
CROSSING_THRESHOLD_LEG = {"gtex_protein_coding": "onebrain", "baron_protein_coding": "oneislets"}
CORR_SUM_POWER = {"corr_sum": 1, "corr_sq_sum": 2}
CROSSING_CATEGORIES = ["above0.5_to_below", "below0.5_to_above", "above0.8_to_below", "below0.8_to_above"]


def base_gct(source):
    return f"{GENE_SPEC}/{source}/{STUB[source]}_protein_coding.gct.gz"


def leg_gct_path(source, leg, agg_suffix=""):
    return f"{GENE_SPEC}/{source}/{leg}{agg_suffix}/{STUB[source]}_{leg}.gct.gz"


def leg_tau_path(source, leg, agg_suffix=""):
    return f"{GENE_SPEC}/{source}/{leg}{agg_suffix}/weighted_tau.csv"


def leg_members(source, leg):
    kind, val = SOURCES[source]["legs"][leg]
    if kind == "prefix":
        return [t for t in TISSUES if t.split("_")[0] == val]
    return val  # kind == "group": explicit member list


LEGS_RE = r"onebrain|oneislets|simplified_0\.1"
AGGS = ["mean", "max"]
AGG_SUFFIX_RE = r"_mean|_max"

wildcard_constraints:
    source = r"gtex_protein_coding|baron_protein_coding",
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
        expand(f"{GENE_SPEC}/{{source}}/figures/tissue_clusters_dendrogram_simplified_0.1.png", source=SOURCES),
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
         for source, cfg in SOURCES.items() for leg in cfg["legs"]]

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
        f"{GENE_SPEC}/{{source}}/figures/tissue_clusters_dendrogram_simplified_0.1.png"
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
        f"{GENE_SPEC}/{{source}}/simplified_0.1{{agg_suffix}}/{{stub}}_simplified_0.1.gct.gz"
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

# ---- corr_sum / corr_sq_sum: closed-form weighting from row sums of the tissue x
# tissue correlation matrix (see gls_tau.corr_sum_weights) -- no matrix inversion, no
# regularization search, unlike diag_reg_inv_corr. corr_sum uses raw r, corr_sq_sum
# uses r^2 (CORR_SUM_POWER); {corr_method} wildcard picks between them rather than a
# separate rule pair, since the two only differ by the power passed to the shared
# compute/plot scripts. Agg/leg-independent, max-agg reference only for now, same as
# diag_reg_inv_corr/margin. ----

rule compute_corr_sum_weighted_tau:
    input:
        lambda wc: base_gct(wc.source)
    output:
        tau=f"{GENE_SPEC}/{{source}}/{{corr_method}}_weights/weighted_tau.csv",
        weights=f"{GENE_SPEC}/{{source}}/{{corr_method}}_weights/weights.csv"
    wildcard_constraints:
        corr_method=r"corr_sum|corr_sq_sum"
    params:
        power=lambda wc: CORR_SUM_POWER[wc.corr_method]
    shell:
        "python3 scripts/compute_corr_sum_weighted_tau.py {input} {params.power} {output.tau} {output.weights}"

rule plot_corr_sum_weights:
    input:
        weights=f"{GENE_SPEC}/{{source}}/{{corr_method}}_weights/weights.csv",
        clusters=lambda wc: CLUSTERS_CSV[wc.source]
    output:
        f"{GENE_SPEC}/{{source}}/{{corr_method}}_weights/figures/tissue_weight_{{corr_method}}.png"
    wildcard_constraints:
        corr_method=r"corr_sum|corr_sq_sum"
    shell:
        "python3 scripts/plot_corr_sum_weights.py {input.weights} {input.clusters} {output} {wildcards.corr_method}"

rule plot_tau_corr_sum_scenarios:
    input:
        corr_sum_tau=f"{GENE_SPEC}/{{source}}/{{corr_method}}_weights/weighted_tau.csv",
        ref_tau=lambda wc: leg_tau_path(wc.source, wc.leg, "_max")
    output:
        f"{GENE_SPEC}/{{source}}/{{corr_method}}_weights/figures/tau_{{corr_method}}_scenarios_scatter_{{leg}}.png"
    wildcard_constraints:
        corr_method=r"corr_sum|corr_sq_sum"
    shell:
        "python3 scripts/plot_tau_scenarios_scatter.py {input.corr_sum_tau} {input.ref_tau} {output} {wildcards.leg} "
        "tau_weighted_{wildcards.corr_method} {wildcards.corr_method}-weighted"
