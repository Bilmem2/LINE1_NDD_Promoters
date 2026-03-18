#!/usr/bin/env python3
"""
05_statistics.py
Performs all statistical tests reported in the manuscript:
  - Mann-Whitney U test for promoter LINE-1 occupancy (each NDD set vs housekeeping)
  - Rank-biserial r effect size with approximate 95% CIs
  - Spearman correlation (intronic length vs LINE-1 density)
  - Length-matched control analysis for intronic confounding

Usage: python scripts/05_statistics.py
Run from the repository root.
"""

import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

RESULTS = Path("results")

# ── Helper functions ──────────────────────────────────────────────────────────

def rank_biserial_r(u_stat, n1, n2):
    """Rank-biserial correlation from Mann-Whitney U statistic."""
    return 1 - (2 * u_stat) / (n1 * n2)

def rb_ci(r, n1, n2, z=1.96):
    """Approximate 95% CI for rank-biserial r."""
    se = np.sqrt((n1 + n2 + 1) / (3 * n1 * n2))
    return (r - z * se, r + z * se)

# ── Load data ─────────────────────────────────────────────────────────────────
print("Loading gene set promoter counts...")
hk = pd.read_csv(RESULTS / "promoter_counts_housekeeping.tsv", sep="\t")

gene_sets = {
    "HPO Seizure":         "hpo_seizure",
    "NDD Tier 1":          "sfari_tier1",
    "HPO ADHD":            "hpo_adhd",
    "Syndromic NDD":       "sfari_syndromic",
    "HPO Seizure ∩ ADHD":  "seizure_adhd",
}

# ── Promoter occupancy analysis ───────────────────────────────────────────────
print("\n=== Promoter LINE-1 Occupancy Analysis ===")
print(f"Housekeeping: n={len(hk)}, occupancy={hk['line1_count'].gt(0).mean()*100:.1f}%\n")

rows = []
for label, fname in gene_sets.items():
    ndd = pd.read_csv(RESULTS / f"promoter_counts_{fname}.tsv", sep="\t")
    u, p = stats.mannwhitneyu(
        ndd["line1_count"], hk["line1_count"], alternative="two-sided"
    )
    r = rank_biserial_r(u, len(ndd), len(hk))
    ci = rb_ci(r, len(ndd), len(hk))
    occ = ndd["line1_count"].gt(0).mean() * 100
    rows.append({
        "Gene Set": label,
        "n": len(ndd),
        "LINE-1 (%)": round(occ, 1),
        "p-value": round(p, 4),
        "r": round(r, 3),
        "CI_lower": round(ci[0], 3),
        "CI_upper": round(ci[1], 3),
    })
    sig = "**" if p < 0.01 else ("*" if p < 0.05 else "n.s.")
    print(f"{label:<28} n={len(ndd):>5}  L1={occ:.1f}%  p={p:.4f}{sig}  r={r:.3f} [{ci[0]:.2f},{ci[1]:.2f}]")

results_df = pd.DataFrame(rows)
results_df.to_csv(RESULTS / "promoter_statistics.tsv", sep="\t", index=False)
print(f"\nSaved: {RESULTS}/promoter_statistics.tsv")

# ── Intronic confounding analysis ─────────────────────────────────────────────
print("\n=== Intronic LINE-1 Confounding Analysis ===")

def gene_level_summary(intron_df):
    """Aggregate intron-level counts to gene level."""
    intron_df["length_bp"] = intron_df["end"] - intron_df["start"]
    gene_summary = intron_df.groupby("gene").agg(
        total_length_kb=("length_bp", lambda x: x.sum() / 1000),
        total_line1=("line1_count", "sum"),
        n_introns=("line1_count", "count"),
    ).reset_index()
    gene_summary["line1_per_kb"] = (
        gene_summary["total_line1"] / gene_summary["total_length_kb"]
    )
    return gene_summary

ndd_intron = pd.read_csv(RESULTS / "intron_counts_sfari_tier1.tsv", sep="\t")
hk_intron  = pd.read_csv(RESULTS / "intron_counts_housekeeping.tsv", sep="\t")

ndd_gene = gene_level_summary(ndd_intron)
hk_gene  = gene_level_summary(hk_intron)

print(f"NDD Tier1 median intronic length: {ndd_gene['total_length_kb'].median():.1f} kb")
print(f"Housekeeping median intronic length: {hk_gene['total_length_kb'].median():.1f} kb")

# Raw comparison
u_raw, p_raw = stats.mannwhitneyu(
    ndd_gene["line1_per_kb"], hk_gene["line1_per_kb"], alternative="two-sided"
)
r_raw = rank_biserial_r(u_raw, len(ndd_gene), len(hk_gene))
print(f"\nRaw intronic: p={p_raw:.2e}, r={r_raw:.3f}")

# Spearman correlation (NDD genes)
rho, p_corr = stats.spearmanr(
    ndd_gene["total_length_kb"], ndd_gene["line1_per_kb"]
)
print(f"Spearman (NDD length vs density): rho={rho:.2f}, p={p_corr:.2e}")

# Length-matched analysis
print("\nRunning length-matched control analysis...")
hk_available = hk_gene.copy().reset_index(drop=True)
matched_pairs = []
used_hk = set()

for _, ndd_row in ndd_gene.iterrows():
    target = ndd_row["total_length_kb"]
    lower, upper = target * 0.5, target * 1.5
    candidates = hk_available[
        (hk_available["total_length_kb"] >= lower) &
        (hk_available["total_length_kb"] <= upper) &
        (~hk_available.index.isin(used_hk))
    ]
    if len(candidates) > 0:
        best_idx = (candidates["total_length_kb"] - target).abs().idxmin()
        matched_pairs.append({
            "ndd_gene": ndd_row["gene"],
            "ndd_l1_per_kb": ndd_row["line1_per_kb"],
            "hk_gene": candidates.loc[best_idx, "gene"],
            "hk_l1_per_kb": candidates.loc[best_idx, "line1_per_kb"],
        })
        used_hk.add(best_idx)

matched_df = pd.DataFrame(matched_pairs)
print(f"Matched pairs: {len(matched_df)}")

u_m, p_m = stats.mannwhitneyu(
    matched_df["ndd_l1_per_kb"],
    matched_df["hk_l1_per_kb"],
    alternative="two-sided"
)
r_m = rank_biserial_r(u_m, len(matched_df), len(matched_df))
print(f"After matching: p={p_m:.3f}, r={r_m:.3f}")

matched_df.to_csv(RESULTS / "intron_matched_pairs.tsv", sep="\t", index=False)
print(f"Saved: {RESULTS}/intron_matched_pairs.tsv")

print("\nStatistical analysis complete.")
