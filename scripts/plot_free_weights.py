"""Horizontal bar chart of the free_weights per-tissue weight (winner of N DE chains,
see optimize_free_weights.py / combine_free_weights_chains.py) for every GTEx tissue,
sorted ascending. Brain_*/Kidney_*/Liver* tissues highlighted separately. Vertical lines
mark the [0.01, 1] optimization bounds (see optimize_free_weights.py)."""
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_WEIGHTS = sys.argv[1]
OUT_PNG = sys.argv[2]

LO_BOUND, HI_BOUND = 0.01, 1.0
CATEGORIES = {"Brain": "tab:red", "Kidney": "tab:green", "Liver": "tab:orange"}

df = pd.read_csv(IN_WEIGHTS).sort_values("weight", ascending=True)

category = pd.Series("other", index=df.index)
for cat in CATEGORIES:
    category[df["tissue"].str.contains(cat)] = cat
colors = category.map({**CATEGORIES, "other": "tab:blue"})

fig, ax = plt.subplots(figsize=(8, 16))
ax.barh(df["tissue"], df["weight"], color=colors)
ax.axvline(1 / len(df), ls="--", color="gray", lw=1, label=f"uniform (1/{len(df)})")
ax.axvline(LO_BOUND, ls=":", color="black", lw=1, label=f"bounds [{LO_BOUND}, {HI_BOUND}]")
ax.axvline(HI_BOUND, ls=":", color="black", lw=1)
ax.set_xlabel("free weight")
ax.set_ylabel("tissue")
ax.tick_params(axis="y", labelsize=6)
ax.legend(handles=[
    plt.Rectangle((0, 0), 1, 1, color="tab:red", label="Brain_*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:green", label="Kidney_*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:orange", label="Liver*"),
    plt.Rectangle((0, 0), 1, 1, color="tab:blue", label="other"),
    ax.lines[0],
    ax.lines[1],
], loc="lower right", fontsize=8)
ax.set_title(f"free_weights: optimized per-tissue weight (n={len(df)} tissues, winner chain)")
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300)
plt.close(fig)
print(f"tissues={len(df)}", file=sys.stderr)
