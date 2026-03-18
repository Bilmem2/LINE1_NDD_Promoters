#!/usr/bin/env bash
# 01_download_data.sh
# Downloads all publicly available raw data sources for the LINE-1 NDD analysis.
# Run from the repository root: bash scripts/01_download_data.sh

set -euo pipefail

mkdir -p data/raw
cd data/raw

echo "=== Downloading UCSC RepeatMasker (hg38) ==="
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz
echo "Decompressing..."
gunzip -k rmsk.txt.gz

echo ""
echo "=== Downloading GENCODE v47 basic annotation GTF ==="
wget -nc https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gtf.gz
echo "Decompressing..."
gunzip -k gencode.v47.basic.annotation.gtf.gz

echo ""
echo "=== Downloading HPO gene lists ==="
wget -nc "https://hpo.jax.org/api/hpo/gene?disease_id=HP:0001250" \
     -O genes_for_HP_0001250.txt \
     --header="Accept: text/plain" || \
wget -nc "https://hpo.jax.org/data/genes-diseases-associations.tsv" \
     -O hpo_gene_associations.tsv

wget -nc "https://hpo.jax.org/api/hpo/gene?disease_id=HP:0007018" \
     -O genes_for_HP_0007018.txt \
     --header="Accept: text/plain"

echo ""
echo "=== Downloading HRT Atlas housekeeping genes ==="
wget -nc "https://housekeeping.unicamp.br/download.php?filename=housekeeping_data/Housekeeping_GenesHuman.csv" \
     -O Housekeeping_GenesHuman.csv

echo ""
echo "=== SFARI Gene: Manual download required ==="
echo "Please visit: https://gene.sfari.org/tools/"
echo "Download the gene CSV and save it to: data/raw/SFARI-Gene_genes_<release>.csv"

echo ""
echo "All automated downloads complete."
