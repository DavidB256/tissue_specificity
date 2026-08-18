"""plot_tau_jackknife_variance_range.py, but tau_range is normalized by each gene's
canonical all-tissue (no tissue dropped) tau value (TAU_COL in weighted_tau.csv) before
plotting -- tau_var is left as-is."""
import sys
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use("Agg")

IN_CSV = sys.argv[1]
IN_TAU = sys.argv[2]
TAU_COL = sys.argv[3]
OUT_PNG = sys.argv[4]

df = pd.read_csv(IN_CSV)[["gene", "tau_var", "tau_range"]]
tau = pd.read_csv(IN_TAU)[["gene", TAU_COL]]
df = df.merge(tau, on="gene").dropna()
df = df[df[TAU_COL] != 0]
df["tau_range_norm"] = df["tau_range"] / df[TAU_COL]

g = sns.jointplot(data=df, x="tau_var", y="tau_range_norm", kind="scatter",
                   joint_kws={"s": 4, "alpha": 0.15}, marginal_kws={"bins": 60})
g.ax_joint.set_xlabel("tau_var (jackknife)")
g.ax_joint.set_ylabel(f"tau_range (jackknife) / {TAU_COL} (all-tissue)")
g.figure.suptitle(f"jackknife tau variance vs normalized range, n={len(df)} genes", y=1.02)
g.savefig(OUT_PNG, dpi=300)
