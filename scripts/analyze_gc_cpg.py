#!/usr/bin/env python3
"""
analyze_gc_cpg.py
Computes GC content and CpG O/E ratio for each promoter window,
then compares these metrics across NDD and housekeeping gene sets.

Usage: python scripts/analyze_gc_cpg.py
Run from the repository root (~/line1_ndd_analysis).
"""

import re
import numpy as np
import pandas as pd
from scipy import stats
from pathlib import Path

FASTA   = Path("data/processed/promoters_2kb.fa")
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

# ── 1. Parse FASTA → per-gene GC% and CpG O/E ───────────────────────────────

def cpg_oe(seq):
    """CpG Observed/Expected = (CpG count * length) / (C count * G count)."""
    seq = seq.upper()
    n   = len(seq)
    if n < 10:
        return np.nan
    c   = seq.count("C")
    g   = seq.count("G")
    cpg = sum(1 for i in range(n - 1) if seq[i] == "C" and seq[i+1] == "G")
    if c == 0 or g == 0:
        return 0.0
    return (cpg * n) / (c * g)

def gc_pct(seq):
    seq = seq.upper()
    n = len(seq)
    if n == 0:
        return np.nan
    return (seq.count("G") + seq.count("C")) / n * 100

print("Parsing FASTA...")
records = {}
current_name = None
current_seq  = []

with open(FASTA) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_name:
                seq = "".join(current_seq)
                records[current_name] = seq
            # bedtools -name uses gene name as header
            current_name = line[1:].split("::")[0]  # gene name before coordinates
            current_seq  = []
        else:
            current_seq.append(line)
    if current_name:
        records[current_name] = "".join(current_seq)

print(f"  Parsed {len(records)} promoter sequences")

rows = []
for gene, seq in records.items():
    rows.append({
        "gene":   gene,
        "gc_pct": gc_pct(seq),
        "cpg_oe": cpg_oe(seq),
    })

seq_df = pd.DataFrame(rows)
seq_df.to_csv(RESULTS / "promoter_gc_cpg_all.tsv", sep="\t", index=False)
print(f"  Saved: results/promoter_gc_cpg_all.tsv")

# ── 2. Load gene sets from existing TSV files ────────────────────────────────

print("\nLoading gene sets from TSV files...")

def genes_from_tsv(path):
    """Extract gene names from our result TSV files."""
    df = pd.read_csv(path, sep="\t")
    return set(df["gene"].dropna())

hk_genes     = genes_from_tsv(RESULTS / "promoter_counts_housekeeping.tsv")
sfari_tier1  = genes_from_tsv(RESULTS / "promoter_counts_sfari_tier1.tsv")
hpo_seizure  = genes_from_tsv(RESULTS / "promoter_counts_hpo_seizure.tsv")
hpo_adhd     = genes_from_tsv(RESULTS / "promoter_counts_hpo_adhd.tsv")
seizure_adhd = genes_from_tsv(RESULTS / "promoter_counts_seizure_adhd.tsv")

syndromic_tsv = RESULTS / "promoter_counts_sfari_syndromic.tsv"
if syndromic_tsv.exists():
    sfari_syndromic = genes_from_tsv(syndromic_tsv)
else:
    sfari_syndromic = set()

gene_sets = {
    "Housekeeping":        hk_genes,
    "HPO Seizure":         hpo_seizure,
    "NDD Tier 1":          sfari_tier1,
    "HPO ADHD":            hpo_adhd,
    "HPO Seizure ∩ ADHD":  seizure_adhd,
}
if sfari_syndromic:
    gene_sets["Syndromic NDD"] = sfari_syndromic

for label, genes in gene_sets.items():
    print(f"  {label}: {len(genes)} genes")

# ── 3. Compare GC% and CpG O/E across gene sets ──────────────────────────────

print("\n=== GC Content and CpG O/E Across Gene Sets ===\n")
print(f"{'Gene Set':<28} {'n':>5}  {'GC% median':>10}  {'CpG O/E median':>14}  {'p (vs HK, GC)':>14}  {'p (vs HK, CpG)':>15}")
print("-" * 95)

hk_df = seq_df[seq_df["gene"].isin(hk_genes)].dropna()

summary_rows = []
for label, genes in gene_sets.items():
    sub = seq_df[seq_df["gene"].isin(genes)].dropna()
    
    if label == "Housekeeping":
        p_gc  = np.nan
        p_cpg = np.nan
    else:
        _, p_gc  = stats.mannwhitneyu(sub["gc_pct"],  hk_df["gc_pct"],  alternative="two-sided")
        _, p_cpg = stats.mannwhitneyu(sub["cpg_oe"],  hk_df["cpg_oe"],  alternative="two-sided")

    med_gc  = sub["gc_pct"].median()
    med_cpg = sub["cpg_oe"].median()

    p_gc_str  = f"{p_gc:.4f}"  if not np.isnan(p_gc)  else "—"
    p_cpg_str = f"{p_cpg:.4f}" if not np.isnan(p_cpg) else "—"

    print(f"{label:<28} {len(sub):>5}  {med_gc:>10.1f}  {med_cpg:>14.3f}  {p_gc_str:>14}  {p_cpg_str:>15}")

    summary_rows.append({
        "Gene Set":        label,
        "n":               len(sub),
        "GC_median":       round(med_gc, 2),
        "CpG_OE_median":   round(med_cpg, 3),
        "p_GC_vs_HK":      round(p_gc, 4)  if not np.isnan(p_gc)  else np.nan,
        "p_CpG_vs_HK":     round(p_cpg, 4) if not np.isnan(p_cpg) else np.nan,
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(RESULTS / "gc_cpg_statistics.tsv", sep="\t", index=False)
print(f"\nSaved: results/gc_cpg_statistics.tsv")

# ── 4. Interpretation hint ───────────────────────────────────────────────────
print("\n=== Interpretation ===")
ndd = seq_df[seq_df["gene"].isin(sfari_tier1)].dropna()
_, p_gc_ndd  = stats.mannwhitneyu(ndd["gc_pct"], hk_df["gc_pct"],  alternative="two-sided")
_, p_cpg_ndd = stats.mannwhitneyu(ndd["cpg_oe"], hk_df["cpg_oe"],  alternative="two-sided")

if p_gc_ndd < 0.05 or p_cpg_ndd < 0.05:
    print("⚠  Significant sequence composition difference detected.")
    print("   GC or CpG O/E differs between NDD Tier1 and housekeeping promoters.")
    print("   Consider adding this as a covariate or noting it as a potential confounder.")
else:
    print("✓  No significant GC or CpG O/E difference between NDD Tier1 and housekeeping.")
    print("   LINE-1 depletion signal is unlikely to be explained by sequence composition alone.")
