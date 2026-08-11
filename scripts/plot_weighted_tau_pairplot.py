"""Seaborn pairplot (scatterplot matrix) of the three tau variants: unweighted,
GLS-weighted, and ridge-regularized GLS-weighted."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

COLS = ["tau_unweighted", "tau_weighted_gls", "tau_weighted_gls_ridge"]

df = pd.read_csv(IN_CSV)[COLS].dropna()

g = sns.pairplot(df, plot_kws={"s": 5, "alpha": 0.15})
g.figure.suptitle(f"tau variants, n={len(df)} genes", y=1.02)
g.savefig(OUT_PNG, dpi=300)
