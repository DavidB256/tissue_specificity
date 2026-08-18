"""Scatter of tau_var vs tau_range (per-gene, across leave-one-tissue-out jackknife
samples, see compute_tau_jackknife_pertissue.py) with marginal histograms."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")

IN_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_CSV)[["tau_var", "tau_range"]].dropna()

g = sns.jointplot(data=df, x="tau_var", y="tau_range", kind="scatter",
                   joint_kws={"s": 4, "alpha": 0.15}, marginal_kws={"bins": 60})
g.ax_joint.set_xlabel("tau_var (jackknife)")
g.ax_joint.set_ylabel("tau_range (jackknife)")
g.figure.suptitle(f"jackknife tau variance vs range, n={len(df)} genes", y=1.02)
g.savefig(OUT_PNG, dpi=300)
