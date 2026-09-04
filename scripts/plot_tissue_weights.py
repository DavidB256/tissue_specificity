"""Horizontal bar chart of ridge-regularized GLS tissue weight (weight_gls_ridge) for
every GTEx tissue, sorted descending. Brain_* tissues highlighted separately."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
OUT_PNG = sys.argv[2]

df = pd.read_csv(IN_WEIGHTS).sort_values("weight_gls_ridge", ascending=True)
is_brain = df["tissue"].str.contains("Brain")
colors = is_brain.map({True: "tab:red", False: "tab:blue"})

fig, ax = plt.subplots(figsize=(8, 16))
ax.barh(df["tissue"], df["weight_gls_ridge"], color=colors)
ax.set_xlabel("weight_gls_ridge")
ax.set_ylabel("tissue")
ax.tick_params(axis="y", labelsize=6)
ax.legend(handles=[
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="Brain_*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="other"),
], loc="lower right", fontsize=8)
ax.set_title(f"Ridge-regularized GLS tissue weight (n={len(df)} tissues, "
             f"lambda={df['ridge_lambda'].iloc[0]:.4g})")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"tissues={len(df)}", file=sys.stderr)
