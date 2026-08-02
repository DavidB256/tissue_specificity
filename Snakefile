DATA = "/ccb/salz4-4/dbass13/expr_data"
GCT = f"{DATA}/GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz"

rule all:
    input:
        f"{DATA}/figures/gmm_grid.png",
        f"{DATA}/figures/gmm_bic_vs_k.png",
        f"{DATA}/figures/expgauss_grid.png",
        f"{DATA}/figures/expgauss_vs_gmm.png",

rule preprocess:
    input:
        GCT
    output:
        f"{DATA}/processed/log10_tpm_long.parquet"
    shell:
        "python3 scripts/preprocess.py {input} {output}"

rule fit_gmm:
    input:
        f"{DATA}/processed/log10_tpm_long.parquet"
    output:
        f"{DATA}/results/gmm_fits.csv"
    shell:
        "python3 scripts/fit_gmm.py {input} {output}"

rule plot_gmm:
    input:
        parquet=f"{DATA}/processed/log10_tpm_long.parquet",
        fits=f"{DATA}/results/gmm_fits.csv"
    output:
        grid=f"{DATA}/figures/gmm_grid.png",
        bic=f"{DATA}/figures/gmm_bic_vs_k.png"
    shell:
        "python3 scripts/plot_gmm.py {input.parquet} {input.fits} {output.grid} {output.bic}"

rule fit_expgauss:
    input:
        f"{DATA}/processed/log10_tpm_long.parquet"
    output:
        f"{DATA}/results/expgauss_fits.csv"
    shell:
        "python3 scripts/fit_expgauss.py {input} {output}"

rule plot_expgauss:
    input:
        parquet=f"{DATA}/processed/log10_tpm_long.parquet",
        fits=f"{DATA}/results/expgauss_fits.csv",
        gmm_fits=f"{DATA}/results/gmm_fits.csv"
    output:
        grid=f"{DATA}/figures/expgauss_grid.png",
        compare=f"{DATA}/figures/expgauss_vs_gmm.png"
    shell:
        "python3 scripts/plot_expgauss.py {input.parquet} {input.fits} {input.gmm_fits} {output.grid} {output.compare}"
