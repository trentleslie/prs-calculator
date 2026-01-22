# Polygenic Risk Score (PRS) Calculator

Calculate personalized polygenic risk scores from whole genome sequencing data using PGS Catalog scores.

## Overview

This project provides tools for:
1. Downloading PGS Catalog scoring files
2. Extracting relevant variants from VCF data
3. Calculating weighted polygenic risk scores
4. Interpreting scores in context of population distributions

## Directory Structure

```
PRS/
├── README.md           # This file
├── pyproject.toml      # Python dependencies
├── src/
│   └── pgs_calculator.py   # Main calculator script
├── scores/             # Downloaded PGS Catalog files
├── results/            # Calculated risk scores
├── docs/
│   ├── PGSC_CALC_QUICKSTART.md   # pgsc_calc single VCF guide
│   ├── PRS_METHODOLOGY_COMPARISON.md  # Method comparison + validation
│   ├── PGS_SETUP_GUIDE.md        # Detailed setup instructions
│   └── screenshots/              # Visual documentation
└── reference/          # Symlink to hg38 reference genome
```

## Quick Start

### Using pgsc_calc (Recommended)

For ancestry-normalized PRS with GRCh38 data, see the **[pgsc_calc Quickstart Guide](docs/PGSC_CALC_QUICKSTART.md)** - complete step-by-step instructions for running on a single VCF file.

### Using Custom Script (Educational/Debugging)

```bash
# Install dependencies
cd Analysis/PRS
uv sync

# Calculate PRS for a specific trait
uv run python src/pgs_calculator.py \
    --vcf ../../Raw_Data/VCF/*.vcf.gz \
    --pgs-file scores/PGS004237_score.txt \
    --output results/

# List available PGS scores
uv run python src/pgs_calculator.py --list-available
```

**Note**: Always validate PGS scores before clinical interpretation. See [PGS Quality Validation](docs/PRS_METHODOLOGY_COMPARISON.md#pgs-quality-validation-findings) for techniques.

## PGS Catalog

The PGS Catalog (https://www.pgscatalog.org/) provides standardized polygenic scores for various traits.

### Current Scores (All LD-Aware)

| PGS ID | Trait | Method | Variants | Citation |
|--------|-------|--------|----------|----------|
| PGS002308 | Type 2 Diabetes | PRS-CSx | 1,259,754 | Ge T et al., Genome Medicine (2022) |
| PGS004034 | Alzheimer's Disease | LDpred2-auto | 1,046,908 | Monti R et al., AJHG (2024) |
| PGS000027 | Body Mass Index | LDpred | 2,100,302 | Khera AV et al., Cell (2019) |
| PGS004237 | Coronary Artery Disease | LDpred | 1,146,511 | Manikpurage HD et al. (2021) |
| PGS005236 | Coronary Artery Disease | AnnoPred | 2,994,055 | Hu J et al., PLoS Comp Bio (2025) |
| PGS001355 | Coronary Artery Disease | AnnoPred | 2,994,055 | Ye Y et al., Circ Genom Prec Med (2021) |

### LD-Aware Method Selection (Critical)

**Always select PGS scores that properly account for linkage disequilibrium (LD)**. Without LD correction, correlated variants are double-counted, producing inaccurate scores.

**Recommended methods** (in order of preference):
1. **LDpred / LDpred2** - Bayesian with LD reference panel (gold standard)
2. **PRS-CS / PRS-CSx** - Bayesian continuous shrinkage
3. **lassosum / SBayesR** - Penalized regression with LD matrix

**Methods to avoid** for accurate estimation:
- Simple sum of GWAS betas (no LD correction)
- P+T without clumping
- "Genome-wide significant variants only" (minimal LD handling)

Download scores from: https://www.pgscatalog.org/downloads/

## Dependencies

- `requests` - PGS Catalog API access
- `pandas` - Data manipulation
- `pysam` / `cyvcf2` - VCF file parsing (optional, for performance)

## Reference Genome

This project uses the hg38 (GRCh38) reference genome, matching the Nucleus sequencing pipeline. The `reference/` directory contains a symlink to the shared reference at `~/genomics_resources/reference/hg38.fa`.

## Results Interpretation

PRS values are reported as:
- **Raw score**: Sum of weighted effect alleles
- **Percentile**: Position relative to population distribution (when available)
- **Z-score**: Standard deviations from population mean

Higher scores generally indicate higher genetic risk, but interpretation requires trait-specific context.

## Ancestry-Normalized Results (January 2026)

Using pgsc_calc v2.0.0 with HGDP+1kGP reference panel (~3,200 samples), normalized against European (EUR) population.

### PRS Summary Table

| Trait | PGS ID | Raw Score | Z-Score | Percentile | Interpretation |
|-------|--------|-----------|---------|------------|----------------|
| **Type 2 Diabetes** | PGS002308 | 0.436 | **+2.76** | **99.7%** | Very High Genetic Risk |
| Alzheimer's Disease | PGS004034 | 1.678 | +0.49 | 64.0% | Average Risk |
| Body Mass Index | PGS000027 | 17.704 | +0.06 | 51.6% | Average |
| Coronary Artery Disease | PGS004237 | -0.258 | -1.47 | 7.2% | Below Average Risk |

### Ancestry Analysis

- **Most Similar Population:** EUR (European)
- **Confidence:** Low (mixed ancestry pattern)
- **RF Probabilities:** EUR 45%, MID 32%, AMR 15%, AFR 3%, CSA 3%, EAS 2%

### Key Findings

1. **Type 2 Diabetes Risk (99.7th percentile):** Exceptionally high genetic predisposition. This aligns with known diagnosis and reinforces the importance of strict glycemic control and lifestyle optimization.

2. **Coronary Artery Disease (7.2nd percentile):** Below average genetic risk. However, T2D significantly increases cardiovascular risk through non-genetic mechanisms - maintain cardiovascular prevention strategies.

3. **BMI (51.6th percentile):** Average genetic predisposition for BMI. Weight management challenges are more likely driven by metabolic factors (diabetes, medications) than genetics.

4. **Alzheimer's Disease (64th percentile):** Moderately elevated. Combined with APOE e3/e4 status, suggests importance of neuroprotective strategies (exercise, metabolic health, cognitive engagement).

### Files

- `results/prs_normalized_ancestry.csv` - Summary table
- `results/pgsc_calc/trent_pgs.txt.gz` - Full scores with all reference samples
- `results/pgsc_calc/trent_popsimilarity.txt.gz` - Detailed ancestry analysis
- `results/pgsc_calc/report.html` - Interactive HTML report
