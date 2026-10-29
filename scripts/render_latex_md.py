"""Render each $$...$$ block of a markdown file as a stacked equation image (matplotlib
mathtext, no TeX install needed). Usage: render_latex_md.py IN_MD OUT_PNG [DPI]"""
import re
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IN_MD, OUT_PNG = sys.argv[1], sys.argv[2]
DPI = int(sys.argv[3]) if len(sys.argv) > 3 else 600
FONTSIZE = 28
ROW_HEIGHT = 1.1  # inches per equation

# mathtext lacks \argmax and \big*; sizes via \left/\right already used elsewhere
SUBS = [(r"\argmax", r"\operatorname{argmax}"), (r"\bigg", ""), (r"\Big", ""), (r"\big", "")]

blocks = [b.strip() for b in re.findall(r"\$\$(.*?)\$\$", open(IN_MD).read(), re.S)]
fig = plt.figure(figsize=(13, ROW_HEIGHT * len(blocks)))
for i, b in enumerate(blocks):
    for old, new in SUBS:
        b = b.replace(old, new)
    fig.text(0.5, 1 - (i + 0.5) / len(blocks), f"${' '.join(b.split())}$",
             ha="center", va="center", fontsize=FONTSIZE)
fig.savefig(OUT_PNG, dpi=DPI, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"{len(blocks)} equations -> {OUT_PNG}", file=sys.stderr)
