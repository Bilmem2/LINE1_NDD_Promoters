#!/usr/bin/env python3
"""
04_gene_sets.py
Curates NDD and housekeeping gene sets, removes overlaps,
and filters the genome-wide LINE-1 count tables by gene set.

Usage: python scripts/04_gene_sets.py
Run from the repository root.
"""

import pandas as pd
from pathlib import Path

# ── Paths ────────────────────────────────────────────────────────────────────
RAW = Path("data/raw")
DATA = Path("data")
RESULTS = Path("results")
RESULTS.mkdir(exist_ok=True)

# ── Load genome-wide LINE-1 counts ───────────────────────────────────────────
print("Loading genome-wide promoter LINE-1 counts...")
cols = ["chr", "start", "end", "gene", "score", "strand", "line1_count"]
promoter_counts = pd.read_csv(
    RESULTS / "promoter_line1_counts_all.bed",
    sep="\t", header=None, names=cols
)

intron_counts = pd.read_csv(
    RESULTS / "intron_line1_counts_all.bed",
    sep="\t", header=None,
    names=["chr", "start", "end", "gene", "score", "strand", "line1_count"]
)

# ── Load gene sets ────────────────────────────────────────────────────────────
print("Loading gene sets...")

# SFARI Gene — adjust filename to match downloaded file
sfari_files = sorted(RAW.glob("SFARI-Gene_genes_*.csv"))
if not sfari_files:
    raise FileNotFoundError(
        "SFARI Gene CSV not found in data/raw/. "
        "Please download from https://gene.sfari.org/tools/"
    )
sfari = pd.read_csv(sfari_files[-1])
sfari_tier1 = set(sfari[sfari["gene-score"].isin([1, 2])]["gene-symbol"].dropna())
sfari_tier2 = set(sfari[sfari["gene-score"].isin([1, 2, 3])]["gene-symbol"].dropna())
sfari_syndromic = set(sfari[sfari["syndromic"] == 1]["gene-symbol"].dropna())
print(f"  SFARI Tier1: {len(sfari_tier1)} genes")
print(f"  SFARI Syndromic: {len(sfari_syndromic)} genes")

# HPO gene lists
hpo_seizure = pd.read_csv(RAW / "genes_for_HP_0001250.txt", sep=r"\s+")
hpo_seizure_genes = set(hpo_seizure["name"].dropna())

hpo_adhd = pd.read_csv(RAW / "genes_for_HP_0007018.txt", sep=r"\s+")
hpo_adhd_genes = set(hpo_adhd["name"].dropna())

seizure_adhd = hpo_seizure_genes & hpo_adhd_genes
print(f"  HPO Seizure: {len(hpo_seizure_genes)} genes")
print(f"  HPO ADHD: {len(hpo_adhd_genes)} genes")
print(f"  Seizure ∩ ADHD: {len(seizure_adhd)} genes")

# HRT Atlas housekeeping genes
hk_raw = pd.read_csv(RAW / "Housekeeping_GenesHuman.csv", sep=";")
hk_genes_raw = set(hk_raw["Gene.name"].dropna())
# Remove genes present in SFARI Tier2 to ensure independence
hk_genes = hk_genes_raw - sfari_tier2
print(f"  Housekeeping (after overlap removal): {len(hk_genes)} genes")

# ── Filter promoter counts by gene set ───────────────────────────────────────
gene_sets = {
    "housekeeping": hk_genes,
    "sfari_tier1": sfari_tier1,
    "sfari_syndromic": sfari_syndromic,
    "hpo_seizure": hpo_seizure_genes,
    "hpo_adhd": hpo_adhd_genes,
    "seizure_adhd": seizure_adhd,
}

print("\nFiltering promoter counts by gene set...")
for name, genes in gene_sets.items():
    subset = promoter_counts[promoter_counts["gene"].isin(genes)].copy()
    out_path = RESULTS / f"promoter_counts_{name}.tsv"
    subset.to_csv(out_path, sep="\t", index=False)
    print(f"  {name}: {len(subset)} promoters → {out_path}")

# ── Filter intron counts for NDD Tier1 and housekeeping ──────────────────────
print("\nFiltering intron counts...")
for name in ["sfari_tier1", "housekeeping"]:
    genes = gene_sets[name]
    subset = intron_counts[intron_counts["gene"].isin(genes)].copy()
    out_path = RESULTS / f"intron_counts_{name}.tsv"
    subset.to_csv(out_path, sep="\t", index=False)
    print(f"  {name}: {len(subset)} intervals → {out_path}")

print("\nGene set curation complete.")
