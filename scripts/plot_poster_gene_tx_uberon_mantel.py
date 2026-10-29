"""Poster bar fig, styled like plot_gene_vs_tx_bars.py (same colors/hatches/
legend/rotated-anchor x labels/font sizes), showing compute_gene_tx_uberon_
mantel.py's Mantel r per scenario (3 sources x TPM<1/TPM<10 floor) with optional bootstrap 95% CI error bars and
significance stars (*** p<0.001, ** p<0.01, * p<0.05) above each bar. dpi=600."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

IN_MANTEL_CSV = sys.argv[1]
OUT_PNG = sys.argv[2]
EXCLUDE = sys.argv[3].split(",") if len(sys.argv) > 3 else []  # optional: sources to drop
IN_BOOT_CSV = sys.argv[4] if len(sys.argv) > 4 else None  # optional: bootstrap 95% CIs (ci_lo/ci_hi per scenario)

SOURCE_ORDER = ["gtex_protein_coding", "gtex_v11_tx_protein_coding", "chess_tx_protein_coding"]
SOURCE_LABEL = {
    "gtex_protein_coding": "GENCODE (gene)",
    "gtex_v11_tx_protein_coding": "GENCODE (tx)",
    "chess_tx_protein_coding": "CHESS (tx)",
}
KIND_LABEL = {"gene_centric": "Gene-level", "tx_centric": "Tx-level"}
KIND_COLOR = {"gene_centric": "#C44E52", "tx_centric": "#4C72B0"}
FLOOR_LABEL = {1.0: "TPM<1 filt.", 10.0: "TPM<10 filt."}
FLOOR_HATCH = {1.0: "", 10.0: "//"}


def short_name(source, floor):
    return f"{SOURCE_LABEL[source]}\n{FLOOR_LABEL[floor]}"


df = pd.read_csv(IN_MANTEL_CSV)
df = df[~df["source"].isin(EXCLUDE)]
df["order_rank"] = df["source"].map({s: i for i, s in enumerate(SOURCE_ORDER)})
df = df.sort_values(["order_rank", "floor_threshold"]).reset_index(drop=True)
df["stars"] = df["stars"].fillna("")

if IN_BOOT_CSV:
    df = df.merge(pd.read_csv(IN_BOOT_CSV)[["scenario", "ci_lo", "ci_hi"]], on="scenario")
else:
    df["ci_lo"] = df["ci_hi"] = df["mantel_r"]

labels = [short_name(s, f) for s, f in zip(df["source"], df["floor_threshold"])]
x = np.arange(len(df))
colors = df["kind"].map(KIND_COLOR)
hatches = df["floor_threshold"].map(FLOOR_HATCH)

fig, ax = plt.subplots(figsize=(8 * len(df) / 6, 6.5))

bars = ax.bar(x, df["mantel_r"], color=colors, edgecolor="black")
for bar, hatch in zip(bars, hatches):
    bar.set_hatch(hatch)
if IN_BOOT_CSV:
    ax.errorbar(x, df["mantel_r"], yerr=[df["mantel_r"] - df["ci_lo"], df["ci_hi"] - df["mantel_r"]],
                fmt="none", ecolor="black", elinewidth=2, capsize=6, capthick=2, zorder=3)
for xi, r, hi, s in zip(x, df["mantel_r"], df["ci_hi"], df["stars"]):
    label = f"{r:.3f} ({s})" if s else f"{r:.3f}"
    ax.text(xi, hi + 0.006, label, ha="center", fontsize=12)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=12, rotation=30, ha="right", rotation_mode="anchor")
ax.set_ylabel("agreement with UBERON (Mantel r of pairwise distances)", fontsize=12)
ax.tick_params(axis="y", labelsize=12)
bottom, top = ax.get_ylim()
ax.set_ylim(bottom, top + 0.1)

legend_elems = [Patch(facecolor=c, edgecolor="black", label=KIND_LABEL[k]) for k, c in KIND_COLOR.items()]
legend_elems += [Patch(facecolor="white", edgecolor="black", hatch=h, label=FLOOR_LABEL[f])
                  for f, h in FLOOR_HATCH.items()]
ax.legend(handles=legend_elems, fontsize=9, loc="upper left")

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=600)
plt.close(fig)

print(f"-> {OUT_PNG}", file=sys.stderr)
