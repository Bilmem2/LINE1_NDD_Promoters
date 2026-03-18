#!/usr/bin/env bash
# 03_bedtools_intersect.sh
# Runs BEDTools intersection to count LINE-1 elements
# overlapping promoter and intronic regions for each gene set.
# Run from the repository root: bash scripts/03_bedtools_intersect.sh

set -euo pipefail

LINE1="data/LINE1_hg38.bed"
PROMOTERS="data/promoters_2kb_hg38.bed"
INTRONS="data/introns_hg38_sorted.bed"
OUTDIR="results"

mkdir -p "$OUTDIR"

echo "=== Counting LINE-1 elements per promoter (genome-wide) ==="
bedtools intersect -c \
    -a "$PROMOTERS" \
    -b "$LINE1" \
    -sorted \
> "$OUTDIR/promoter_line1_counts_all.bed"
echo "Done: $OUTDIR/promoter_line1_counts_all.bed"

echo ""
echo "=== Counting LINE-1 elements per intronic interval (genome-wide) ==="
bedtools intersect -c \
    -a "$INTRONS" \
    -b "$LINE1" \
    -sorted \
> "$OUTDIR/intron_line1_counts_all.bed"
echo "Done: $OUTDIR/intron_line1_counts_all.bed"

echo ""
echo "BEDTools intersection complete."
echo "Run scripts/04_gene_sets.py to filter by gene set."
