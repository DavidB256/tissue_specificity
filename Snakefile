import sys
sys.path.insert(0, "scripts")
import pandas as pd
from model_registry import MODELS

DATA = "/ccb/salz4-4/dbass13/expr_data"
GCT = f"{DATA}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"
MODEL_NAMES = list(MODELS.keys())
TISSUES = list(pd.read_csv(GCT, sep="\t", skiprows=2, nrows=0).columns[2:])

rule all:
    input:
        f"{DATA}/results/model_comparison.csv",
        f"{DATA}/methods.md",
        expand(f"{DATA}/figures/{{model}}_grid.png", model=MODEL_NAMES)

rule preprocess:
    input:
        GCT
    output:
        f"{DATA}/processed/log10_tpm_long.parquet"
    shell:
        "python3 scripts/preprocess.py {input} {output}"

rule fit_one:
    input:
        f"{DATA}/processed/log10_tpm_long.parquet"
    output:
        f"{DATA}/results/by_model_tissue/{{model}}/{{tissue}}.csv"
    shell:
        "python3 scripts/fit_one.py {wildcards.model} {wildcards.tissue} {input} {output}"

rule combine_model:
    input:
        lambda wc: expand(f"{DATA}/results/by_model_tissue/{wc.model}/{{tissue}}.csv", tissue=TISSUES)
    output:
        f"{DATA}/results/{{model}}_fits.csv"
    run:
        dfs = [pd.read_csv(f) for f in input]
        pd.concat(dfs, ignore_index=True).to_csv(output[0], index=False)

rule combine_comparison:
    input:
        expand(f"{DATA}/results/{{model}}_fits.csv", model=MODEL_NAMES)
    output:
        f"{DATA}/results/model_comparison.csv"
    shell:
        "python3 scripts/combine_comparison.py {DATA}/results {output}"

rule plot_model_grid:
    input:
        parquet=f"{DATA}/processed/log10_tpm_long.parquet",
        fits=f"{DATA}/results/{{model}}_fits.csv"
    output:
        f"{DATA}/figures/{{model}}_grid.png"
    shell:
        "python3 scripts/plot_model_grid.py {wildcards.model} {input.parquet} {input.fits} {output}"

rule write_methods:
    output:
        f"{DATA}/methods.md"
    shell:
        "python3 scripts/write_methods.py {output}"
