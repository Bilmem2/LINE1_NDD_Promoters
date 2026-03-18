#!/usr/bin/env bash
# 02_prepare_regions.sh
# Extracts promoter windows and intronic intervals from GENCODE v47.
# Filters LINE-1 elements from UCSC RepeatMasker annotation.
# Run from the repository root: bash scripts/02_prepare_regions.sh

set -euo pipefail

GTF="data/raw/gencode.v47.basic.annotation.gtf"
RMSK="data/raw/rmsk.txt"
OUTDIR="data"

echo "=== Extracting LINE-1 elements from RepeatMasker ==="
awk 'NR > 1 && $12 == "LINE" && $13 == "L1" {
    print $6 "\t" $7-1 "\t" $8 "\t" $13 "\t0\t" $10
}' "$RMSK" \
| sort -k1,1 -k2,2n \
> "$OUTDIR/LINE1_hg38.bed"
echo "LINE-1 elements: $(wc -l < "$OUTDIR/LINE1_hg38.bed")"

echo ""
echo "=== Extracting ±2 kb promoter windows (strand-aware) ==="
awk '
$3 == "gene" && $0 ~ /gene_type "protein_coding"/ {
    # Parse chromosome and strand
    chr = $1
    strand = $7
    if (strand == "+") {
        tss = $4 - 1   # convert to 0-based
    } else {
        tss = $5        # end coordinate for - strand
    }
    start = tss - 2000
    end   = tss + 2000
    if (start < 0) start = 0
    # Extract gene_name
    match($0, /gene_name "([^"]+)"/, arr)
    gene = arr[1]
    print chr "\t" start "\t" end "\t" gene "\t0\t" strand
}' "$GTF" \
| sort -k1,1 -k2,2n \
> "$OUTDIR/promoters_2kb_hg38.bed"
echo "Promoter windows: $(wc -l < "$OUTDIR/promoters_2kb_hg38.bed")"

echo ""
echo "=== Extracting intronic intervals ==="
python3 scripts/extract_introns.py \
    --gtf "$GTF" \
    --out "$OUTDIR/introns_hg38_sorted.bed"
echo "Intronic intervals: $(wc -l < "$OUTDIR/introns_hg38_sorted.bed")"

echo ""
echo "Region preparation complete."
