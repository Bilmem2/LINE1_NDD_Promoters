# LINE-1 Depletion at Promoters of Neurodevelopmental Disorder Genes: A Genome-Wide Analysis

> **Note:** This repository serves as the computational supplement for the Brief Report: *"LINE-1 Depletion at Promoters of Neurodevelopmental Disorder Genes: A Genome-Wide Analysis"*.

**Author:** Can Sevilmiş  
**Affiliation:** Department of Molecular Biology and Genetics, Bahçeşehir University, Istanbul, Turkey  
**ORCID:** [0000-0002-9180-1924](https://orcid.org/0000-0002-9180-1924)

---

## Overview

This repository contains all scripts, data processing pipelines, and the manuscript source for the study described above. A genome-wide computational analysis was performed comparing LINE-1 retrotransposon density in promoter regions (±2 kb from TSS) across five independently curated neurodevelopmental disorder (NDD) gene sets against a housekeeping gene control, using GRCh38/hg38 annotations.

**Key Findings:** * NDD gene promoters are consistently depleted of LINE-1 elements relative to housekeeping gene promoters (n=1,982). 
* The magnitude of depletion tracks with phenotypic specificity (Housekeeping > HPO Seizure > NDD Tier 1 > Syndromic NDD > HPO ADHD > HPO Seizure ∩ ADHD).
* **Robustness:** The depletion signal is maintained after controlling for intronic length differences (length-matched analysis: n=678 pairs, p=0.036). Furthermore, no significant differences in GC content (p=0.3289) or CpG density (p=0.9665) were detected between NDD Tier 1 and housekeeping promoters, indicating that the depletion is independent of primary sequence composition biases.

---

## Repository Structure

```
LINE1_NDD_Promoters/
├── manuscript/
│   ├── main.tex          # LaTeX source
│   └── references.bib    # BibTeX reference file
├── scripts/
│   ├── 01_download_data.sh       # Download all raw data sources
│   ├── 02_prepare_regions.sh     # Extract promoter and intron BED files
│   ├── extract_introns.py        # Helper script called by 02 to calculate introns
│   ├── 03_bedtools_intersect.sh  # LINE-1 intersection analysis
│   ├── 04_gene_sets.py           # Gene set curation and overlap removal
│   ├── 05_statistics.py          # Statistical tests and effect sizes
│   ├── 06_figures.py             # Generate all figures
│   └── analyze_gc_cpg.py         # GC content and CpG O/E robustness analysis
├── data/
│   └── README.md            # Data sources and download instructions
├── figures/
│   ├── final_panel_v2.png           # Main figure (promoter occupancy + effect sizes)
│   ├── final_panel_v2.svg           # Vector version of the main figure
│   └── line1_silencing_diagram.png  # Conceptual model diagram
├── results/
│   └── README.md            # Description of result files
├── environment/
│   └── environment.yml      # Conda environment specification
└── LICENSE                  # MIT License
```

---

## Reproducibility

### 1. Set up the environment

```bash
conda env create -f environment/environment.yml
conda activate line1-ndd
```

### 2. Download raw data

```bash
bash scripts/01_download_data.sh
```
*Note: Some external databases (like HRT Atlas or HPO) may experience SSL certificate issues or link changes. If automated downloads fail, please check `data/README.md` for manual download instructions.*

This will download:
- UCSC RepeatMasker hg38 annotation (`rmsk.txt.gz`)
- GENCODE v47 comprehensive annotation GTF
- SFARI Gene list (manual download required — see `data/README.md`)
- HPO gene lists for HP:0001250 and HP:0007018
- HRT Atlas housekeeping gene list

### 3. Prepare genomic regions

```bash
bash scripts/02_prepare_regions.sh
```

Outputs: `data/promoters_2kb_hg38.bed`, `data/introns_hg38_sorted.bed`, `data/LINE1_hg38.bed`

### 4. Run BEDTools intersection

```bash
bash scripts/03_bedtools_intersect.sh
```

### 5. Curate gene sets and run statistics

```bash
python scripts/04_gene_sets.py
python scripts/05_statistics.py
```

### 6. Generate figures

```bash
python scripts/06_figures.py
```

### 7. Run sequence composition robustness analysis (GC/CpG)

To run the nucleotide composition analysis, the human reference genome (hg38) is required.

```bash
# Download hg38 reference genome
wget -nc https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz -O data/raw/hg38.fa.gz
gunzip -k data/raw/hg38.fa.gz

# Extract FASTA sequences and run analysis
bedtools getfasta -fi data/raw/hg38.fa -bed data/promoters_2kb_hg38.bed -nameOnly > data/processed/promoters_2kb.fa
python scripts/analyze_gc_cpg.py
```

---

## Data Sources

| Resource | Version | URL |
|---|---|---|
| UCSC RepeatMasker | hg38 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ |
| GENCODE | v47 | https://www.gencodegenes.org |
| SFARI Gene | March 2026 | https://gene.sfari.org |
| HPO | 2021 release | https://hpo.jax.org |
| HRT Atlas | v1.0 | https://housekeeping.unicamp.br |
| UCSC Reference Genome | hg38 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/ |

---

## Citation

> Sevilmiş C. (2026). Promoter-Specific Depletion of LINE-1 in Neurodevelopmental Disorder Genes.

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
