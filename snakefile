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
GENE_SPEC_GTEX = f"{GENE_SPEC}/gtex"

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

# Salmon quantification-error simulation: does deviation from the exp_skewgauss2 fit
# detect quantification errors? "control" (real reads, unmodified) vs "miscall_k"
# (k% of read positions randomly changed to a different base), k in 1/5 (%).
# Reads pooled from the 10 largest Kidney CRAMs (read-only source, never modified) --
# single-CRAM data gave a degenerate single-peak-near-0 TPM distribution, too little
# data. guided quantification = mapping-based salmon quant against a GENCODE v50
# BASIC transcriptome, decoy-aware index (genome as decoy) built from the same
# reference used to decode the CRAMs -- standard current-best-practice salmon
# indexing; only --recoverOrphans and -p are used as explicit salmon quant flags,
# per instruction (-p bumped to 20 to use more of this machine's cores).
SALMON_SIM = f"{DATA}/salmon_error_sim"
KIDNEY_CRAMS = [  # 10 largest Kidney CRAMs by size, hardcoded (not wildcarded, per instruction)
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR810007.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1416516.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1325483.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1469746.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR821356.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1380931.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1486080.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1432650.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1445835.cram",
    "/ccb/salz8-1/avaraby/GTEx_2019_hisat2/Kidney/SRR1443092.cram",
]
REF_GENOME = "/ccb/salz4-4/dbass13/recursive_introns_data/GRCh38/hg38.p12.fa"
GENCODE_V50_BASIC_URL = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.basic.annotation.gtf.gz"
MISCALL_RATES = ["1", "5"]
SALMON_CONDITIONS = ["control"] + [f"miscall_{k}" for k in MISCALL_RATES]

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
MANE_GTF_SRC = "/ccb/salz4-4/dbass13/recursive_introns_data/annotations/MANE.v1.5.gtf"
GENCODE_V47_GTF = "/ccb/salz4-4/dbass13/recursive_introns_data/annotations/gencode.v47.annotation.gtf"
GENE_SPEC_BARON_SIMPLIFIED = f"{GENE_SPEC_BARON}/simplified_0.1"
BARON_SIMPLIFIED_GCT = f"{GENE_SPEC_BARON_SIMPLIFIED}/celltype_tpm_simplified_0.1.gct.gz"
GENE_SPEC_BARON_FREE = f"{GENE_SPEC_BARON}/free_weights/simplified_0.1"
N_BARON_FREE_WEIGHTS_CHAINS = 10  # 14-dim DE converges fast/consistently (single chain: rho=0.989 in 50s) -- far fewer chains needed than the 68-dim GTEx free_weights leg

wildcard_constraints:
    species = "macaque|marmoset",
    k = r"1|5",
    tissue = r"[^/]+"

rule all:
    input:
        f"{MODEL_FITS}/model_comparison.csv",
        "methods.md",
        expand(f"{DIST_FIGS}/{{model}}_grid.png", model=MODEL_NAMES),
        expand(f"{PCA_FIGS}/{{model}}_pca.png", model=MODEL_NAMES),
        f"{MODEL_FITS}/pca_compressibility_comparison.csv",
        f"{MODEL_FITS}/pca_compressibility_comparison_allgenes.csv",
        expand(f"{DATA}/{{species}}/raw/gene_tpm_combined.gct.gz", species=NHP_LIST),
        expand(f"{DATA}/{{species}}/model_fits/model_comparison.csv", species=NHP_LIST),
        expand(f"{DATA}/{{species}}/distribution_fits/{{model}}_grid.png", species=NHP_LIST, model=MODEL_NAMES),
        f"{SALMON_SIM}/figures/exp_skewgauss2_conditions.png",
        f"{GENE_SPEC_GTEX}/figures/tau_expbeta_scatter.png",
        f"{GENE_SPEC_GTEX}/weighted_tau.csv",
        f"{GENE_SPEC_GTEX}/figures/weighted_tau_pairplot.png",
        f"{GENE_SPEC_GTEX}/figures/weighted_tau_pairplot_squared.png",
        f"{GENE_SPEC_GTEX}/figures/tissue_weight_gls_ridge_squared.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/default/weighted_tau.csv",
        f"{GENE_SPEC_GTEX}/cluster_weights/onebrain/weighted_tau.csv",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tau_cluster_onebrain_vs_default.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tau_cluster_onebrain_vs_default_bkl.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tissue_clusters_pca.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tissue_clusters_dendrogram_default.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tissue_clusters_dendrogram_onebrain.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/param_search_best.csv",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tau_cluster_onebrain_vs_default_naive.png",
        f"{GENE_SPEC_GTEX}/cluster_weights/figures/tau_cluster_onebrain_vs_default_optimized.png",
        f"{GENE_SPEC_GTEX}/free_weights/onebrain/weighted_tau.csv",
        f"{GENE_SPEC_GTEX}/free_weights/simplified_0.1/weighted_tau.csv",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tissue_weight_free_onebrain.png",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tissue_weight_free_simplified_0.1.png",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tau_free_vs_default_onebrain.png",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tau_free_vs_default_simplified_0.1.png",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tau_delta_scatter.png",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tau_delta_scatter_abs.png",
        f"{GENE_SPEC_GTEX}/free_weights/figures/tau_delta_delta_violin.png",
        f"{GENE_SPEC_GTEX}/free_weights/top_delta_delta_tau.csv",
        f"{GENE_SPEC_GTEX}/figures/tissue_weight_gls_ridge.png",
        f"{GENE_SPEC_GTEX}/tau_stability_spearman.csv",
        f"{GENE_SPEC_GTEX}/tau_stability_loto.csv",
        f"{GENE_SPEC_GTEX}/figures/neural_tau_impact.png",
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_clustermap.png",
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_heatmap_only.png",
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_clustermap_squared.png",
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_clustermap_rawtpm.png",
        f"{GENE_SPEC_GTEX}/figures/tau_neural_downsample.png",
        f"{GENE_SPEC_BARON}/celltype_tpm.csv",
        f"{GENE_SPEC_BARON}/figures/tissue_correlation_clustermap.png",
        f"{GENE_SPEC_BARON}/figures/tissue_clusters_dendrogram_simplified_0.1.png",
        f"{GENE_SPEC_BARON}/figures/tissue_weight_free_simplified_0.1_by_cluster.png",
        f"{GENE_SPEC_BARON}/figures/tau_scenarios_scatter_simplified_0.1.png"

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

rule download_gencode_v50_basic:
    output:
        f"{SALMON_SIM}/ref/gencode.v50.basic.annotation.gtf.gz"
    params:
        url=GENCODE_V50_BASIC_URL
    shell:
        "curl -sL -o {output} '{params.url}'"

rule extract_reads_from_cram:
    input:
        KIDNEY_CRAMS
    output:
        r1=f"{SALMON_SIM}/reads/control_R1.fastq.gz",
        r2=f"{SALMON_SIM}/reads/control_R2.fastq.gz"
    params:
        genome=REF_GENOME,
        tmpdir=f"{SALMON_SIM}/reads/tmp_percram"
    threads: 16
    shell:
        """
        mkdir -p {params.tmpdir}
        rm -f {output.r1} {output.r2}
        for cram in {input}; do
            name=$(basename "$cram" .cram)
            samtools collate -f -@ 16 --reference {params.genome} -O "$cram" \
                | samtools fastq -@ 16 -1 {params.tmpdir}/${{name}}_R1.fastq.gz -2 {params.tmpdir}/${{name}}_R2.fastq.gz -0 /dev/null -s /dev/null -F 0x900 -c 1 -
            cat {params.tmpdir}/${{name}}_R1.fastq.gz >> {output.r1}
            cat {params.tmpdir}/${{name}}_R2.fastq.gz >> {output.r2}
            rm -f {params.tmpdir}/${{name}}_R1.fastq.gz {params.tmpdir}/${{name}}_R2.fastq.gz
        done
        rmdir {params.tmpdir} 2>/dev/null || true
        """

rule simulate_miscalls:
    input:
        r1=f"{SALMON_SIM}/reads/control_R1.fastq.gz",
        r2=f"{SALMON_SIM}/reads/control_R2.fastq.gz"
    output:
        r1=f"{SALMON_SIM}/reads/miscall_{{k}}_R1.fastq",
        r2=f"{SALMON_SIM}/reads/miscall_{{k}}_R2.fastq"
    shell:
        "python3 scripts/simulate_miscalls.py {input.r1} {input.r2} {wildcards.k} {output.r1} {output.r2} 0"

rule decompress_gtf:
    input:
        f"{SALMON_SIM}/ref/gencode.v50.basic.annotation.gtf.gz"
    output:
        f"{SALMON_SIM}/ref/gencode.v50.basic.annotation.gtf"
    shell:
        "gunzip -k -c {input} > {output}"

rule extract_transcripts_fa:
    input:
        gtf=f"{SALMON_SIM}/ref/gencode.v50.basic.annotation.gtf",
        genome=REF_GENOME
    output:
        f"{SALMON_SIM}/ref/transcripts.fa"
    shell:
        "/ccb/sw/bin/gffread -w {output} -g {input.genome} {input.gtf}"

rule build_decoys:
    input:
        REF_GENOME
    output:
        f"{SALMON_SIM}/ref/decoys.txt"
    shell:
        "grep '^>' {input} | cut -d ' ' -f1 | sed 's/>//' > {output}"

rule build_gentrome:
    input:
        transcripts=f"{SALMON_SIM}/ref/transcripts.fa",
        genome=REF_GENOME
    output:
        f"{SALMON_SIM}/ref/gentrome.fa"
    shell:
        "cat {input.transcripts} {input.genome} > {output}"

rule build_salmon_index:
    input:
        gentrome=f"{SALMON_SIM}/ref/gentrome.fa",
        decoys=f"{SALMON_SIM}/ref/decoys.txt"
    output:
        directory(f"{SALMON_SIM}/ref/salmon_index")
    threads: 8
    shell:
        "salmon index -t {input.gentrome} -d {input.decoys} -i {output} -p 8"

rule salmon_quant_control:
    input:
        r1=f"{SALMON_SIM}/reads/control_R1.fastq.gz",
        r2=f"{SALMON_SIM}/reads/control_R2.fastq.gz",
        index=f"{SALMON_SIM}/ref/salmon_index"
    output:
        f"{SALMON_SIM}/quant/control/quant.sf"
    params:
        outdir=f"{SALMON_SIM}/quant/control"
    threads: 20
    shell:
        "salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -p 20 --recoverOrphans -o {params.outdir}"

rule salmon_quant_miscall:
    input:
        r1=f"{SALMON_SIM}/reads/miscall_{{k}}_R1.fastq",
        r2=f"{SALMON_SIM}/reads/miscall_{{k}}_R2.fastq",
        index=f"{SALMON_SIM}/ref/salmon_index"
    output:
        f"{SALMON_SIM}/quant/miscall_{{k}}/quant.sf"
    params:
        outdir=f"{SALMON_SIM}/quant/miscall_{{k}}"
    threads: 20
    shell:
        "salmon quant -i {input.index} -l A -1 {input.r1} -2 {input.r2} -p 20 --recoverOrphans -o {params.outdir}"

rule preprocess_salmon_conditions:
    input:
        expand(f"{SALMON_SIM}/quant/{{cond}}/quant.sf", cond=SALMON_CONDITIONS)
    output:
        f"{SALMON_SIM}/processed/salmon_conditions_long.parquet"
    params:
        args=lambda wc, input: " ".join(f"{c} {f}" for c, f in zip(SALMON_CONDITIONS, input))
    shell:
        "python3 scripts/preprocess_salmon_conditions.py {params.args} {output}"

rule fit_exp_skewgauss2_salmon_condition:
    input:
        f"{SALMON_SIM}/processed/salmon_conditions_long.parquet"
    output:
        f"{SALMON_SIM}/model_fits/by_condition/exp_skewgauss2/{{condition}}.csv"
    shell:
        "python3 scripts/fit_one.py exp_skewgauss2 {wildcards.condition} {input} {output}"

rule combine_salmon_condition_fits:
    input:
        expand(f"{SALMON_SIM}/model_fits/by_condition/exp_skewgauss2/{{condition}}.csv", condition=SALMON_CONDITIONS)
    output:
        f"{SALMON_SIM}/model_fits/exp_skewgauss2_fits.csv"
    run:
        dfs = [pd.read_csv(f) for f in input]
        pd.concat(dfs, ignore_index=True).to_csv(output[0], index=False)

rule plot_salmon_conditions:
    input:
        parquet=f"{SALMON_SIM}/processed/salmon_conditions_long.parquet",
        fits=f"{SALMON_SIM}/model_fits/exp_skewgauss2_fits.csv"
    output:
        f"{SALMON_SIM}/figures/exp_skewgauss2_conditions.png"
    shell:
        "python3 scripts/plot_salmon_conditions.py {input.parquet} {input.fits} {output}"

rule compute_tau:
    input:
        GCT
    output:
        f"{GENE_SPEC_GTEX}/tau.csv"
    shell:
        "python3 scripts/compute_tau.py {input} {output}"

rule compute_weighted_tau:
    input:
        GCT
    output:
        tau=f"{GENE_SPEC_GTEX}/weighted_tau.csv",
        weights=f"{GENE_SPEC_GTEX}/weighted_tau_weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau.py {input} {output.tau} {output.weights}"

rule compute_tau_stability:
    input:
        gct=GCT,
        tau=f"{GENE_SPEC_GTEX}/weighted_tau.csv"
    output:
        spearman=f"{GENE_SPEC_GTEX}/tau_stability_spearman.csv",
        loto=f"{GENE_SPEC_GTEX}/tau_stability_loto.csv"
    shell:
        "python3 scripts/compute_tau_stability.py {input.gct} {input.tau} {output.spearman} {output.loto}"

rule identify_neural_genes:
    input:
        GCT
    output:
        f"{GENE_SPEC_GTEX}/neural_genes.csv"
    shell:
        "python3 scripts/identify_neural_genes.py {input} {output}"

rule plot_neural_tau_impact:
    input:
        tau=f"{GENE_SPEC_GTEX}/weighted_tau.csv",
        neural=f"{GENE_SPEC_GTEX}/neural_genes.csv"
    output:
        csv=f"{GENE_SPEC_GTEX}/neural_tau_impact.csv",
        png=f"{GENE_SPEC_GTEX}/figures/neural_tau_impact.png"
    shell:
        "python3 scripts/plot_neural_tau_impact.py {input.tau} {input.neural} {output.csv} {output.png}"

rule plot_weighted_tau_pairplot:
    input:
        f"{GENE_SPEC_GTEX}/weighted_tau.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/weighted_tau_pairplot.png"
    shell:
        "python3 scripts/plot_weighted_tau_pairplot.py {input} {output}"

rule compute_weighted_tau_squared:
    input:
        GCT
    output:
        tau=f"{GENE_SPEC_GTEX}/weighted_tau_squared.csv",
        weights=f"{GENE_SPEC_GTEX}/weighted_tau_squared_weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau_squared.py {input} {output.tau} {output.weights}"

rule plot_weighted_tau_pairplot_squared:
    input:
        f"{GENE_SPEC_GTEX}/weighted_tau_squared.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/weighted_tau_pairplot_squared.png"
    shell:
        "python3 scripts/plot_weighted_tau_pairplot.py {input} {output}"

GENE_SPEC_ONEBRAIN = f"{GENE_SPEC_GTEX}/onebrain"
ONEBRAIN_GCT = f"{GENE_SPEC_ONEBRAIN}/gene_median_tpm_onebrain.gct.gz"

rule build_gct_onetissue:
    input:
        GCT
    output:
        f"{GENE_SPEC_GTEX}/one{{tissue}}/gene_median_tpm_one{{tissue}}.gct.gz"
    params:
        prefix=lambda wildcards: wildcards.tissue.capitalize()
    shell:
        "python3 scripts/build_gct_onetissue.py {input} {output} {params.prefix}"

rule compute_tissue_correlation_onebrain:
    input:
        ONEBRAIN_GCT
    output:
        corr=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation_linkage.csv"
    shell:
        "python3 scripts/compute_tissue_correlation.py {input} {output.corr} {output.linkage}"

rule compute_weighted_tau_onebrain:
    input:
        ONEBRAIN_GCT
    output:
        tau=f"{GENE_SPEC_ONEBRAIN}/weighted_tau.csv",
        weights=f"{GENE_SPEC_ONEBRAIN}/weighted_tau_weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau.py {input} {output.tau} {output.weights}"

GENE_SPEC_CLUSTER_WEIGHTS = f"{GENE_SPEC_GTEX}/cluster_weights"

rule compute_tau_cluster_weighted_default:
    input:
        gct=GCT,
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation.csv"
    output:
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/tissue_clusters.csv",
        tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/weighted_tau.csv"
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} {output.clusters} {output.tau}"

rule compute_tau_cluster_weighted_onebrain:
    input:
        gct=ONEBRAIN_GCT,
        linkage=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation.csv"
    output:
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain/tissue_clusters.csv",
        tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain/weighted_tau.csv"
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} {output.clusters} {output.tau}"

rule plot_tau_cluster_vs_default:
    input:
        default_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/weighted_tau.csv",
        onebrain_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain/weighted_tau.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tau_cluster_onebrain_vs_default.png"
    shell:
        "python3 scripts/plot_tau_cluster_vs_default.py {input.default_tau} {input.onebrain_tau} {output} onebrain"

rule plot_tau_cluster_vs_default_bkl:
    input:
        gct=GCT,
        default_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/weighted_tau.csv",
        onebrain_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain/weighted_tau.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tau_cluster_onebrain_vs_default_bkl.png"
    shell:
        "python3 scripts/plot_tau_cluster_vs_default_bkl.py {input.gct} {input.default_tau} {input.onebrain_tau} {output} onebrain"

rule plot_tissue_clusters_pca:
    input:
        gct=GCT,
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/tissue_clusters.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tissue_clusters_pca.png"
    shell:
        "python3 scripts/plot_tissue_clusters_pca.py {input.gct} {input.clusters} {output}"

rule plot_tissue_clusters_dendrogram_default:
    input:
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv",
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/tissue_clusters.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tissue_clusters_dendrogram_default.png"
    shell:
        "python3 scripts/plot_tissue_clusters_dendrogram.py {input.linkage} {input.clusters} {output} default"

rule plot_tissue_clusters_dendrogram_onebrain:
    input:
        linkage=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation_linkage.csv",
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain/tissue_clusters.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tissue_clusters_dendrogram_onebrain.png"
    shell:
        "python3 scripts/plot_tissue_clusters_dendrogram.py {input.linkage} {input.clusters} {output} onebrain"

rule optimize_cluster_weight_params:
    input:
        default_gct=GCT,
        default_linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv",
        default_corr=f"{GENE_SPEC_GTEX}/tissue_correlation.csv",
        onebrain_gct=ONEBRAIN_GCT,
        onebrain_linkage=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation_linkage.csv",
        onebrain_corr=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation.csv"
    output:
        trajectory=f"{GENE_SPEC_CLUSTER_WEIGHTS}/param_search_trajectory.csv",
        best=f"{GENE_SPEC_CLUSTER_WEIGHTS}/param_search_best.csv"
    shell:
        "python3 scripts/optimize_cluster_weight_params.py {input.default_gct} {input.default_linkage} "
        "{input.default_corr} {input.onebrain_gct} {input.onebrain_linkage} {input.onebrain_corr} "
        "{output.trajectory} {output.best}"

# best params from optimize_cluster_weight_params (param_search_best.csv), pinned as constants
# so the naive-vs-optimized comparison figures don't silently drift if the search is rerun.
DIST_THRESHOLD_OPT = "0.013112"
WEIGHT_CAP_OPT = "0.100032"

rule compute_tau_cluster_weighted_default_optimized:
    input:
        gct=GCT,
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation.csv"
    output:
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default_optimized/tissue_clusters.csv",
        tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default_optimized/weighted_tau.csv"
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} "
        "{output.clusters} {output.tau} " + DIST_THRESHOLD_OPT + " " + WEIGHT_CAP_OPT

rule compute_tau_cluster_weighted_onebrain_optimized:
    input:
        gct=ONEBRAIN_GCT,
        linkage=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_ONEBRAIN}/tissue_correlation.csv"
    output:
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain_optimized/tissue_clusters.csv",
        tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain_optimized/weighted_tau.csv"
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} "
        "{output.clusters} {output.tau} " + DIST_THRESHOLD_OPT + " " + WEIGHT_CAP_OPT

rule plot_tau_cluster_vs_default_naive:
    input:
        default_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/weighted_tau.csv",
        onebrain_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain/weighted_tau.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tau_cluster_onebrain_vs_default_naive.png"
    shell:
        "python3 scripts/plot_tau_cluster_vs_default.py {input.default_tau} {input.onebrain_tau} {output} "
        "'onebrain (t=0.10, c=0.50)'"

rule plot_tau_cluster_vs_default_optimized:
    input:
        default_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default_optimized/weighted_tau.csv",
        onebrain_tau=f"{GENE_SPEC_CLUSTER_WEIGHTS}/onebrain_optimized/weighted_tau.csv"
    output:
        f"{GENE_SPEC_CLUSTER_WEIGHTS}/figures/tau_cluster_onebrain_vs_default_optimized.png"
    shell:
        "python3 scripts/plot_tau_cluster_vs_default.py {input.default_tau} {input.onebrain_tau} {output} "
        "'onebrain (t=0.013, c=0.10)'"

GENE_SPEC_CLUSTER_SIMPLIFIED = f"{GENE_SPEC_GTEX}/simplified_0.1"
CLUSTER_SIMPLIFIED_GCT = f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/gene_median_tpm_simplified_0.1.gct.gz"

rule build_gct_cluster_simplified:
    input:
        gct=GCT,
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/tissue_clusters.csv"
    output:
        CLUSTER_SIMPLIFIED_GCT
    shell:
        "python3 scripts/build_gct_clustered.py {input.gct} {input.clusters} {output}"

rule compute_weighted_tau_cluster_simplified:
    input:
        CLUSTER_SIMPLIFIED_GCT
    output:
        tau=f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau.csv",
        weights=f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau_weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau.py {input} {output.tau} {output.weights}"

# free_weights: one directory per collapse scenario ("leg"), each holding its own
# independently-optimized 68-dim weight vector (see optimize_free_weights.py) fit
# against that scenario's own reduced-resolution reference tau. Figures for both legs
# live side by side in free_weights/figures/.
GENE_SPEC_FREE_WEIGHTS = f"{GENE_SPEC_GTEX}/free_weights"
N_FREE_WEIGHTS_CHAINS = 60

FREE_WEIGHTS_REF_TAU = {
    "onebrain": f"{GENE_SPEC_ONEBRAIN}/weighted_tau.csv",
    "simplified_0.1": f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau.csv",
}

def free_weights_ref_tau(wildcards):
    return FREE_WEIGHTS_REF_TAU[wildcards.leg]

rule optimize_free_weights_chain:
    input:
        default_gct=GCT,
        ref_tau=free_weights_ref_tau
    output:
        weights=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/chains/chain_{{chain}}/weights.csv",
        tau=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/chains/chain_{{chain}}/weighted_tau.csv",
        log=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/chains/chain_{{chain}}/optimization_log.csv",
        best=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/chains/chain_{{chain}}/best.csv"
    wildcard_constraints:
        leg=r"onebrain|simplified_0\.1"
    shell:
        "python3 scripts/optimize_free_weights.py {input.default_gct} {input.ref_tau} "
        "{output.weights} {output.tau} {output.log} {output.best} {wildcards.chain}"

rule combine_free_weights_chains:
    input:
        lambda wildcards: expand(
            f"{GENE_SPEC_FREE_WEIGHTS}/{wildcards.leg}/chains/chain_{{chain}}/best.csv",
            chain=range(N_FREE_WEIGHTS_CHAINS))
    output:
        weights=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/weights.csv",
        tau=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/weighted_tau.csv",
        summary=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/chains_summary.csv"
    wildcard_constraints:
        leg=r"onebrain|simplified_0\.1"
    shell:
        "python3 scripts/combine_free_weights_chains.py {GENE_SPEC_FREE_WEIGHTS}/{wildcards.leg}/chains "
        "{N_FREE_WEIGHTS_CHAINS} {output.weights} {output.tau} {output.summary}"

rule plot_free_weights:
    input:
        f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/weights.csv"
    output:
        f"{GENE_SPEC_FREE_WEIGHTS}/figures/tissue_weight_free_{{leg}}.png"
    wildcard_constraints:
        leg=r"onebrain|simplified_0\.1"
    shell:
        "python3 scripts/plot_free_weights.py {input} {output}"

rule plot_tau_free_vs_default:
    input:
        gct=GCT,
        free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/{{leg}}/weighted_tau.csv",
        ref_tau=free_weights_ref_tau
    output:
        f"{GENE_SPEC_FREE_WEIGHTS}/figures/tau_free_vs_default_{{leg}}.png"
    wildcard_constraints:
        leg=r"onebrain|simplified_0\.1"
    shell:
        "python3 scripts/plot_tau_free_vs_default.py {input.gct} {input.free_tau} {input.ref_tau} {output} {wildcards.leg}"

rule plot_tau_delta_scatter:
    input:
        onebrain_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/onebrain/weighted_tau.csv",
        onebrain_ref_tau=f"{GENE_SPEC_ONEBRAIN}/weighted_tau.csv",
        simplified_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/simplified_0.1/weighted_tau.csv",
        simplified_ref_tau=f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau.csv"
    output:
        f"{GENE_SPEC_FREE_WEIGHTS}/figures/tau_delta_scatter.png"
    shell:
        "python3 scripts/plot_tau_delta_scatter.py {input.onebrain_free_tau} {input.onebrain_ref_tau} "
        "{input.simplified_free_tau} {input.simplified_ref_tau} {output}"

rule plot_tau_delta_scatter_abs:
    input:
        onebrain_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/onebrain/weighted_tau.csv",
        onebrain_ref_tau=f"{GENE_SPEC_ONEBRAIN}/weighted_tau.csv",
        simplified_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/simplified_0.1/weighted_tau.csv",
        simplified_ref_tau=f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau.csv"
    output:
        f"{GENE_SPEC_FREE_WEIGHTS}/figures/tau_delta_scatter_abs.png"
    shell:
        "python3 scripts/plot_tau_delta_scatter_abs.py {input.onebrain_free_tau} {input.onebrain_ref_tau} "
        "{input.simplified_free_tau} {input.simplified_ref_tau} {output}"

rule plot_tau_delta_delta_violin:
    input:
        gct=GCT,
        clusters=f"{GENE_SPEC_CLUSTER_WEIGHTS}/default/tissue_clusters.csv",
        onebrain_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/onebrain/weighted_tau.csv",
        onebrain_ref_tau=f"{GENE_SPEC_ONEBRAIN}/weighted_tau.csv",
        simplified_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/simplified_0.1/weighted_tau.csv",
        simplified_ref_tau=f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau.csv"
    output:
        f"{GENE_SPEC_FREE_WEIGHTS}/figures/tau_delta_delta_violin.png"
    shell:
        "python3 scripts/plot_tau_delta_delta_violin.py {input.gct} {input.clusters} "
        "{input.onebrain_free_tau} {input.onebrain_ref_tau} "
        "{input.simplified_free_tau} {input.simplified_ref_tau} {output}"

rule make_top_delta_delta_tau_table:
    input:
        gct=GCT,
        onebrain_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/onebrain/weighted_tau.csv",
        onebrain_ref_tau=f"{GENE_SPEC_ONEBRAIN}/weighted_tau.csv",
        simplified_free_tau=f"{GENE_SPEC_FREE_WEIGHTS}/simplified_0.1/weighted_tau.csv",
        simplified_ref_tau=f"{GENE_SPEC_CLUSTER_SIMPLIFIED}/weighted_tau.csv"
    output:
        f"{GENE_SPEC_FREE_WEIGHTS}/top_delta_delta_tau.csv"
    shell:
        "python3 scripts/make_top_delta_delta_tau_table.py {input.gct} "
        "{input.onebrain_free_tau} {input.onebrain_ref_tau} "
        "{input.simplified_free_tau} {input.simplified_ref_tau} {output}"

rule plot_tissue_weights_squared:
    input:
        f"{GENE_SPEC_GTEX}/weighted_tau_squared_weights.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tissue_weight_gls_ridge_squared.png"
    shell:
        "python3 scripts/plot_tissue_weights.py {input} {output}"

rule plot_tissue_weights:
    input:
        f"{GENE_SPEC_GTEX}/weighted_tau_weights.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tissue_weight_gls_ridge.png"
    shell:
        "python3 scripts/plot_tissue_weights.py {input} {output}"

rule compute_tissue_correlation:
    input:
        GCT
    output:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv"
    shell:
        "python3 scripts/compute_tissue_correlation.py {input} {output.corr} {output.linkage}"

rule plot_tissue_correlation:
    input:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_clustermap.png"
    shell:
        "python3 scripts/plot_tissue_correlation.py {input.corr} {input.linkage} {output}"

rule plot_tissue_correlation_heatmap_only:
    input:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_linkage.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_heatmap_only.png"
    shell:
        """python3 scripts/plot_tissue_correlation_heatmap_only.py {input.corr} {input.linkage} {output} \
            "Correlation matrix of log10 TPM values in GTEx tissues" """

rule compute_tissue_correlation_squared:
    input:
        f"{GENE_SPEC_GTEX}/tissue_correlation.csv"
    output:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation_squared.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_squared_linkage.csv"
    shell:
        "python3 scripts/compute_tissue_correlation_squared.py {input} {output.corr} {output.linkage}"

rule plot_tissue_correlation_squared:
    input:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation_squared.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_squared_linkage.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_clustermap_squared.png"
    shell:
        """python3 scripts/plot_tissue_correlation.py {input.corr} {input.linkage} {output} \
            "log10 TPM, TPM<1->0, R elementwise squared" """

rule compute_tissue_correlation_rawtpm:
    input:
        GCT
    output:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation_rawtpm.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_rawtpm_linkage.csv"
    shell:
        "python3 scripts/compute_tissue_correlation_rawtpm.py {input} {output.corr} {output.linkage}"

rule plot_tissue_correlation_rawtpm:
    input:
        corr=f"{GENE_SPEC_GTEX}/tissue_correlation_rawtpm.csv",
        linkage=f"{GENE_SPEC_GTEX}/tissue_correlation_rawtpm_linkage.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tissue_correlation_clustermap_rawtpm.png"
    shell:
        """python3 scripts/plot_tissue_correlation.py {input.corr} {input.linkage} {output} \
            "log10 TPM, no floor (raw TPM)" """

rule compute_tau_neural_downsample:
    input:
        gct=GCT
    output:
        f"{GENE_SPEC_GTEX}/tau_neural_downsample.csv"
    shell:
        "python3 scripts/compute_tau_neural_downsample.py {input.gct} {output}"

rule plot_tau_neural_downsample:
    input:
        tau=f"{GENE_SPEC_GTEX}/weighted_tau.csv",
        downsample=f"{GENE_SPEC_GTEX}/tau_neural_downsample.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tau_neural_downsample.png"
    shell:
        "python3 scripts/plot_tau_neural_downsample.py {input.tau} {input.downsample} {output}"

rule compute_exp_beta_ratio:
    input:
        gct=GCT,
        fits=f"{MODEL_FITS}/exp_beta_fits.csv"
    output:
        f"{GENE_SPEC_GTEX}/exp_beta_ratio.csv"
    shell:
        "python3 scripts/compute_exp_beta_ratio.py {input.gct} {input.fits} {output}"

rule compute_beta_posterior:
    input:
        gct=GCT,
        fits=f"{MODEL_FITS}/exp_beta_fits.csv"
    output:
        f"{GENE_SPEC_GTEX}/beta_posterior.csv"
    shell:
        "python3 scripts/compute_beta_posterior.py {input.gct} {input.fits} {output}"

rule plot_tau_exp_beta:
    input:
        tau=f"{GENE_SPEC_GTEX}/tau.csv",
        posterior=f"{GENE_SPEC_GTEX}/beta_posterior.csv"
    output:
        f"{GENE_SPEC_GTEX}/figures/tau_expbeta_scatter.png"
    shell:
        "python3 scripts/plot_tau_exp_beta.py {input.tau} {input.posterior} {output}"


rule link_mane_gtf:
    output:
        f"{DATA}/MANE.v1.5.gtf"
    shell:
        "ln -s {MANE_GTF_SRC} {output}"

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

rule compute_tissue_correlation_baron:
    input:
        f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz"
    output:
        corr=f"{GENE_SPEC_BARON}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC_BARON}/tissue_correlation_linkage.csv"
    shell:
        "python3 scripts/compute_tissue_correlation.py {input} {output.corr} {output.linkage}"

rule plot_tissue_correlation_baron:
    input:
        corr=f"{GENE_SPEC_BARON}/tissue_correlation.csv",
        linkage=f"{GENE_SPEC_BARON}/tissue_correlation_linkage.csv"
    output:
        f"{GENE_SPEC_BARON}/figures/tissue_correlation_clustermap.png"
    shell:
        """python3 scripts/plot_tissue_correlation.py {input.corr} {input.linkage} {output} \
            "log10 TPM, TPM<1->0, Baron pancreas cell types" """

rule compute_tau_cluster_weighted_baron:
    input:
        gct=f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz",
        linkage=f"{GENE_SPEC_BARON}/tissue_correlation_linkage.csv",
        corr=f"{GENE_SPEC_BARON}/tissue_correlation.csv"
    output:
        clusters=f"{GENE_SPEC_BARON}/cluster_weights/tissue_clusters_simplified_0.1.csv",
        tau=f"{GENE_SPEC_BARON}/cluster_weights/weighted_tau_simplified_0.1.csv"
    shell:
        "python3 scripts/compute_tau_cluster_weighted.py {input.gct} {input.linkage} {input.corr} {output.clusters} {output.tau}"

rule plot_tissue_clusters_dendrogram_baron:
    input:
        linkage=f"{GENE_SPEC_BARON}/tissue_correlation_linkage.csv",
        clusters=f"{GENE_SPEC_BARON}/cluster_weights/tissue_clusters_simplified_0.1.csv"
    output:
        f"{GENE_SPEC_BARON}/figures/tissue_clusters_dendrogram_simplified_0.1.png"
    shell:
        "python3 scripts/plot_tissue_clusters_dendrogram.py {input.linkage} {input.clusters} {output} 'Baron pancreas cell types'"

rule build_baron_gct_simplified:
    input:
        gct=f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz",
        clusters=f"{GENE_SPEC_BARON}/cluster_weights/tissue_clusters_simplified_0.1.csv"
    output:
        BARON_SIMPLIFIED_GCT
    shell:
        "python3 scripts/build_gct_clustered.py {input.gct} {input.clusters} {output}"

rule compute_weighted_tau_baron_simplified:
    input:
        BARON_SIMPLIFIED_GCT
    output:
        tau=f"{GENE_SPEC_BARON_SIMPLIFIED}/weighted_tau.csv",
        weights=f"{GENE_SPEC_BARON_SIMPLIFIED}/weighted_tau_weights.csv"
    shell:
        "python3 scripts/compute_weighted_tau.py {input} {output.tau} {output.weights}"

rule optimize_baron_free_weights_chain:
    input:
        default_gct=f"{GENE_SPEC_BARON}/celltype_tpm.gct.gz",
        ref_tau=f"{GENE_SPEC_BARON_SIMPLIFIED}/weighted_tau.csv"
    output:
        weights=f"{GENE_SPEC_BARON_FREE}/chains/chain_{{chain}}/weights.csv",
        tau=f"{GENE_SPEC_BARON_FREE}/chains/chain_{{chain}}/weighted_tau.csv",
        log=f"{GENE_SPEC_BARON_FREE}/chains/chain_{{chain}}/optimization_log.csv",
        best=f"{GENE_SPEC_BARON_FREE}/chains/chain_{{chain}}/best.csv"
    shell:
        "python3 scripts/optimize_free_weights.py {input.default_gct} {input.ref_tau} "
        "{output.weights} {output.tau} {output.log} {output.best} {wildcards.chain}"

rule combine_baron_free_weights_chains:
    input:
        expand(f"{GENE_SPEC_BARON_FREE}/chains/chain_{{chain}}/best.csv", chain=range(N_BARON_FREE_WEIGHTS_CHAINS))
    output:
        weights=f"{GENE_SPEC_BARON_FREE}/weights.csv",
        tau=f"{GENE_SPEC_BARON_FREE}/weighted_tau.csv",
        summary=f"{GENE_SPEC_BARON_FREE}/chains_summary.csv"
    shell:
        "python3 scripts/combine_free_weights_chains.py {GENE_SPEC_BARON_FREE}/chains "
        "{N_BARON_FREE_WEIGHTS_CHAINS} {output.weights} {output.tau} {output.summary}"

rule plot_baron_free_weights_by_cluster:
    input:
        weights=f"{GENE_SPEC_BARON_FREE}/weights.csv",
        clusters=f"{GENE_SPEC_BARON}/cluster_weights/tissue_clusters_simplified_0.1.csv"
    output:
        f"{GENE_SPEC_BARON}/figures/tissue_weight_free_simplified_0.1_by_cluster.png"
    shell:
        "python3 scripts/plot_free_weights_by_cluster.py {input.weights} {input.clusters} {output}"

rule plot_baron_tau_scenarios_scatter:
    input:
        free_tau=f"{GENE_SPEC_BARON_FREE}/weighted_tau.csv",
        ref_tau=f"{GENE_SPEC_BARON_SIMPLIFIED}/weighted_tau.csv"
    output:
        f"{GENE_SPEC_BARON}/figures/tau_scenarios_scatter_simplified_0.1.png"
    shell:
        "python3 scripts/plot_tau_scenarios_scatter.py {input.free_tau} {input.ref_tau} {output} simplified_0.1"

