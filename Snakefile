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
UMAP_FIGS = f"{DATA}/tissue_umap"
UMAP_MODELS = ["exp_skewgauss2"]  # trying UMAP out on this model first

# dGTEx non-human-primate data: gene_tpm.gct.gz per species is per-SAMPLE, not
# per-tissue-median like the human file. The bucket also publishes a pre-split
# tpms-by-tissue/ directory (per-tissue files, samples as columns) -- used here to
# compute the per-gene-per-tissue median directly, rather than reimplementing the
# sample->tissue mapping ourselves.
NHP_SPECIES = {
    "macaque": {
        "url_dir": "macaques",
        "token": "macaque",
        "combined_url": "https://storage.googleapis.com/nhp-dgtex/bulk-gex/v1/rna-seq/macaques/NHP-dGTEx_Macaque_Analysis_2025-07-21_v1_RNASeQCv2.4.3_gene_tpm.gct.gz",
        "tissues": ["adipose_tissue", "adrenal_gland", "blood_vessel", "brain", "colon",
                    "esophagus", "heart", "kidney", "liver", "lung", "lymph_node",
                    "muscle", "nerve", "ovary", "pancreas", "pituitary", "skin",
                    "small_intestine", "spleen", "stomach", "testis", "thymus",
                    "thyroid", "trachea", "uterus", "vagina"],
    },
    "marmoset": {
        "url_dir": "marmosets",
        "token": "marmosets",
        "combined_url": "https://storage.googleapis.com/nhp-dgtex/bulk-gex/v1/rna-seq/marmosets/NHP-dGTEx_Marmoset_Analysis_2025-08-08_v1_RNASeQCv2.4.3_gene_tpm.gct.gz",
        "tissues": ["adipose_tissue", "adrenal_gland", "blood_vessel", "brain", "colon",
                    "esophagus", "heart", "kidney", "liver", "lung", "lymph_node",
                    "muscle", "nerve", "ovary", "pancreas", "pituitary", "skin",
                    "small_intestine", "spleen", "stomach", "testis", "thymus",
                    "thyroid", "trachea", "uterus"],  # no vagina file published
    },
}
NHP_LIST = list(NHP_SPECIES.keys())

wildcard_constraints:
    species = "macaque|marmoset"

rule all:
    input:
        f"{MODEL_FITS}/model_comparison.csv",
        f"{DATA}/methods.md",
        expand(f"{DIST_FIGS}/{{model}}_grid.png", model=MODEL_NAMES),
        expand(f"{PCA_FIGS}/{{model}}_pca.png", model=MODEL_NAMES),
        expand(f"{UMAP_FIGS}/{{model}}_umap.png", model=UMAP_MODELS),
        f"{MODEL_FITS}/pca_compressibility_comparison.csv",
        f"{MODEL_FITS}/pca_compressibility_comparison_allgenes.csv",
        expand(f"{DATA}/{{species}}/raw/gene_tpm_combined.gct.gz", species=NHP_LIST),
        expand(f"{DATA}/{{species}}/model_fits/model_comparison.csv", species=NHP_LIST),
        expand(f"{DATA}/{{species}}/distribution_fits/{{model}}_grid.png", species=NHP_LIST, model=MODEL_NAMES)

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
    run:
        dfs = [pd.read_csv(f) for f in input]
        pd.concat(dfs, ignore_index=True).to_csv(output[0], index=False)

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

rule download_nhp_combined:
    output:
        f"{DATA}/{{species}}/raw/gene_tpm_combined.gct.gz"
    params:
        url=lambda wc: NHP_SPECIES[wc.species]["combined_url"]
    shell:
        "curl -sL -o {output} '{params.url}'"

rule download_nhp_tissue:
    output:
        f"{DATA}/{{species}}/raw/tissue_{{tissue}}.gct.gz"
    params:
        url=lambda wc: (f"https://storage.googleapis.com/nhp-dgtex/bulk-gex/v1/rna-seq/"
                         f"{NHP_SPECIES[wc.species]['url_dir']}/tpms-by-tissue/"
                         f"gene_tpm_nhp-{NHP_SPECIES[wc.species]['token']}_v1_{wc.tissue}.gct.gz")
    shell:
        "curl -sL -o {output} '{params.url}'"

rule build_nhp_median_gct:
    input:
        lambda wc: [f"{DATA}/{wc.species}/raw/tissue_{t}.gct.gz" for t in NHP_SPECIES[wc.species]["tissues"]]
    output:
        f"{DATA}/{{species}}/processed/gene_median_tpm.gct.gz"
    params:
        tissues=lambda wc: ",".join(NHP_SPECIES[wc.species]["tissues"])
    shell:
        "python3 scripts/build_median_gct.py {params.tissues} {input} {output}"

rule nhp_preprocess:
    input:
        f"{DATA}/{{species}}/processed/gene_median_tpm.gct.gz"
    output:
        f"{DATA}/{{species}}/processed/log10_tpm_long.parquet"
    shell:
        "python3 scripts/preprocess.py {input} {output}"

rule nhp_fit_one:
    input:
        f"{DATA}/{{species}}/processed/log10_tpm_long.parquet"
    output:
        f"{DATA}/{{species}}/model_fits/by_model_tissue/{{model}}/{{tissue}}.csv"
    shell:
        "python3 scripts/fit_one.py {wildcards.model} {wildcards.tissue} {input} {output}"

rule nhp_combine_model:
    input:
        lambda wc: expand(f"{DATA}/{wc.species}/model_fits/by_model_tissue/{wc.model}/{{tissue}}.csv",
                           tissue=NHP_SPECIES[wc.species]["tissues"])
    output:
        f"{DATA}/{{species}}/model_fits/{{model}}_fits.csv"
    run:
        dfs = [pd.read_csv(f) for f in input]
        pd.concat(dfs, ignore_index=True).to_csv(output[0], index=False)

rule nhp_combine_comparison:
    input:
        lambda wc: expand(f"{DATA}/{wc.species}/model_fits/{{model}}_fits.csv", model=MODEL_NAMES)
    output:
        f"{DATA}/{{species}}/model_fits/model_comparison.csv"
    params:
        fits_dir=lambda wc: f"{DATA}/{wc.species}/model_fits"
    shell:
        "python3 scripts/combine_comparison.py {params.fits_dir} {output}"

rule nhp_plot_model_grid:
    input:
        parquet=f"{DATA}/{{species}}/processed/log10_tpm_long.parquet",
        fits=f"{DATA}/{{species}}/model_fits/{{model}}_fits.csv"
    output:
        f"{DATA}/{{species}}/distribution_fits/{{model}}_grid.png"
    shell:
        "python3 scripts/plot_model_grid.py {wildcards.model} {input.parquet} {input.fits} {output}"

rule pca_tissue_clustering:
    input:
        f"{MODEL_FITS}/{{model}}_fits.csv"
    output:
        f"{PCA_FIGS}/{{model}}_pca.png"
    shell:
        "python3 scripts/pca_tissue_clustering.py {wildcards.model} {input} {output}"

rule umap_tissue_clustering:
    input:
        f"{MODEL_FITS}/{{model}}_fits.csv"
    output:
        f"{UMAP_FIGS}/{{model}}_umap.png"
    shell:
        "python3 scripts/umap_tissue_clustering.py {wildcards.model} {input} {output}"

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
        f"{DATA}/methods.md"
    shell:
        "python3 scripts/write_methods.py {output}"
