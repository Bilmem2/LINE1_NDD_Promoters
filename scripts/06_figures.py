#!/usr/bin/env python3
"""
06_figures.py
Generates Figure 1 (final_panel_v2): bar chart of LINE-1 promoter occupancy
and effect size forest plot across NDD gene sets vs housekeeping.

Usage: python scripts/06_figures.py
Run from the repository root.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from pathlib import Path

RESULTS = Path("results")
FIGURES = Path("figures")
FIGURES.mkdir(exist_ok=True)

# ── Load statistics table ─────────────────────────────────────────────────────
stats_df = pd.read_csv(RESULTS / "promoter_statistics.tsv", sep="\t")

# Add housekeeping baseline row
hk_row = pd.DataFrame([{
    "Gene Set": "Housekeeping",
    "n": 1982,
    "LINE-1 (%)": 31.1,
    "p-value": np.nan,
    "r": 0.0,
    "CI_lower": 0.0,
    "CI_upper": 0.0,
}])
plot_df = pd.concat([hk_row, stats_df], ignore_index=True)

# Order for x-axis (from least to most depleted among NDD)
order = [
    "Housekeeping",
    "HPO Seizure",
    "NDD Tier 1",
    "HPO ADHD",
    "Syndromic NDD",
    "HPO Seizure ∩ ADHD",
]
plot_df = plot_df.set_index("Gene Set").loc[order].reset_index()

# ── Significance labels ───────────────────────────────────────────────────────
def sig_label(p):
    if pd.isna(p):
        return ""
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."

plot_df["sig"] = plot_df["p-value"].apply(sig_label)

# ── Colors ────────────────────────────────────────────────────────────────────
colors = ["#95a5a6"] + ["#2c7bb6"] * 5  # grey for HK, blue for NDD

# ── Figure layout ─────────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))

# Panel A: Bar chart
ax = axes[0]
bars = ax.bar(
    range(len(plot_df)),
    plot_df["LINE-1 (%)"],
    color=colors,
    edgecolor="white",
    width=0.65,
)
ax.axhline(31.1, color="#7f8c8d", linestyle="--", linewidth=1.2, zorder=0)
ax.set_xticks(range(len(plot_df)))
ax.set_xticklabels(
    [f"{row['Gene Set']}\n(n={int(row['n'])})" for _, row in plot_df.iterrows()],
    fontsize=8.5, rotation=20, ha="right"
)
ax.set_ylabel("Promoters with ≥1 LINE-1 element (%)", fontsize=10)
ax.set_ylim(0, 38)
ax.set_title("A", fontweight="bold", loc="left", fontsize=12)

# Add significance markers
for i, row in plot_df.iterrows():
    if row["sig"]:
        ax.text(i, row["LINE-1 (%)"] + 0.5, row["sig"],
                ha="center", va="bottom", fontsize=9, fontweight="bold")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

# Panel B: Forest plot (effect sizes)
ax2 = axes[1]
ndd_df = plot_df[plot_df["Gene Set"] != "Housekeeping"].reset_index(drop=True)

y_pos = list(range(len(ndd_df)))
labels = [f"{row['Gene Set']} (n={int(row['n'])})" for _, row in ndd_df.iterrows()]

ax2.errorbar(
    ndd_df["r"],
    y_pos,
    xerr=[ndd_df["r"] - ndd_df["CI_lower"], ndd_df["CI_upper"] - ndd_df["r"]],
    fmt="o",
    color="#2c7bb6",
    ecolor="#5b9ec9",
    capsize=4,
    markersize=6,
    linewidth=1.5,
)
ax2.axvline(0, color="#7f8c8d", linestyle="--", linewidth=1)
ax2.set_yticks(y_pos)
ax2.set_yticklabels(labels, fontsize=8.5)
ax2.set_xlabel("Rank-biserial r (95% CI)", fontsize=10)
ax2.set_title("B", fontweight="bold", loc="left", fontsize=12)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

plt.tight_layout(pad=2.0)
out_path = FIGURES / "final_panel_v2.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.savefig(str(out_path).replace(".png", ".svg"), bbox_inches="tight")
print(f"Saved: {out_path}")
print(f"Saved: {str(out_path).replace('.png', '.svg')}")
