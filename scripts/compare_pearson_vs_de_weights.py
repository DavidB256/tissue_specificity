"""Compare free_cluster_weights (DE, Spearman) vs free_cluster_weights_pearson (L-BFGS-B, Pearson r)
per source/leg: Pearson r and Spearman rho of each solution's tau vs the leg's reference tau,
weight agreement (weights scaled to max=1; tau is scale invariant), chain spread, runtimes.
Usage: GENE_SPEC DE_TIMING_DIR OUT_CSV  (DE_TIMING_DIR/<source>_<leg>/runtime.txt = one DE chain wall s)"""
import sys
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

GENE_SPEC, DE_TIMING_DIR, OUT_CSV = sys.argv[1:4]
RUNS = [("gtex_protein_coding", "onebrain"), ("gtex_protein_coding", "hierclust0.1"),
        ("gtex_v11_tx_protein_coding", "onebrain"),
        ("chess_tx_protein_coding", "onebrain"), ("chess_tx_protein_coding", "hierclust0.1"),
        ("baron_protein_coding", "oneislets"), ("baron_protein_coding", "hierclust0.1")]


def metrics(tau_csv, ref):
    t = pd.read_csv(tau_csv).set_index("gene")["tau_free_weighted"]
    d = pd.concat([ref, t], axis=1, join="inner").dropna()
    return pearsonr(d.iloc[:, 0], d.iloc[:, 1])[0], spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0]


rows = []
for source, leg in RUNS:
    base = f"{GENE_SPEC}/{source}"
    ref = pd.read_csv(f"{base}/{leg}_max/weighted_tau.csv").set_index("gene")["tau_unweighted"]
    de, pe = f"{base}/free_cluster_weights_max/{leg}", f"{base}/free_cluster_weights_pearson_max/{leg}"
    r_de, rho_de = metrics(f"{de}/weighted_tau.csv", ref)
    r_pe, rho_pe = metrics(f"{pe}/weighted_tau.csv", ref)
    w_de = pd.read_csv(f"{de}/weights.csv").set_index("tissue")["weight"]
    w_pe = pd.read_csv(f"{pe}/weights.csv").set_index("tissue")["weight"].loc[w_de.index]
    a, b = w_de / w_de.max(), w_pe / w_pe.max()
    cs = pd.read_csv(f"{pe}/chains_summary.csv")
    ds = pd.read_csv(f"{de}/chains_summary.csv")
    de_chain_s = float(open(f"{DE_TIMING_DIR}/{source}_{leg}/runtime.txt").read())
    rows.append(dict(
        source=source, leg=leg, n_tissues=len(w_de),
        r_de=r_de, r_pearson=r_pe, d_r=r_pe - r_de,
        rho_de=rho_de, rho_pearson=rho_pe, d_rho=rho_pe - rho_de,
        w_corr_pearson=pearsonr(a, b)[0], w_corr_spearman=spearmanr(a, b)[0],
        w_max_abs_diff=(a - b).abs().max(), w_mean_abs_diff=(a - b).abs().mean(),
        de_chain_spread_rho=ds["best_rho"].max() - ds["best_rho"].min(),
        pearson_chain_spread_r=cs["best_r"].max() - cs["best_r"].min(),
        de_chain_s=de_chain_s, n_de_chains=len(ds),
        pearson_chain_s_median=cs["runtime_s"].median(), pearson_chain_s_max=cs["runtime_s"].max(),
        n_pearson_chains=len(cs), pearson_starts_per_chain=int(cs["n_starts"].iloc[0]),
        speedup_per_chain=de_chain_s / cs["runtime_s"].median(),
    ))
out = pd.DataFrame(rows)
out.to_csv(OUT_CSV, index=False)
pd.set_option("display.width", 250, "display.max_columns", 50)
print(out.T.to_string(float_format=lambda v: f"{v:.6g}"))
