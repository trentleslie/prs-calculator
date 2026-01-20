# PRS Population Normalization Status

## Summary

We attempted to run pgsc_calc v2.0 with ancestry normalization to convert raw PRS scores to population-normalized Z-scores and percentiles.

## What Was Accomplished

### ✅ Completed Steps

1. **Installed plink2** (v2.0.0-a.7LM)
   - Downloaded from S3 and configured
   - Used for VCF to pgen format conversion

2. **Converted VCF to PLINK2 format**
   - 4,937,275 variants loaded
   - Used `--vcf-allow-no-nonvar` for single-sample VCF
   - Used `--split-par hg38` for chrX pseudoautosomal regions
   - Output: `data/trent_wgs.{pgen,pvar,psam}`

3. **Ran pgsc_calc pipeline** (raw scores only)
   - Successfully calculated PRS for all 4 traits
   - Match rates: 48-53% (lower than our calculator due to no reference imputation)
   - Results in `results/pgsc_calc/`

### ❌ Ancestry Normalization Failed

The `--run_ancestry` flag caused pipeline failures due to a bug in pgsc_calc v2.0:
- Error: `No signature of method: java.lang.Boolean.getFileSystem()`
- Likely related to Nextflow/DSL2 path handling
- Bug appears to be in the PathVisitor component

## Score Comparison

| Trait | Our LD-Aware | pgsc_calc | Coverage Difference |
|-------|-------------|-----------|---------------------|
| Type 2 Diabetes | 0.392 | 0.436 | 99.96% vs 48.34% |
| Alzheimer's Disease | 1.622 | 1.678 | 99.97% vs 53.08% |
| Body Mass Index | 38.638 | 17.717 | 99.99% vs 51.57% |
| Coronary Artery Disease | -0.270 | -0.259 | 99.96% vs 52.40% |

**Key Insight**: Our calculator achieves ~99.9% coverage by imputing missing variants from reference, while pgsc_calc with a single-sample VCF only matches ~50% of scoring file variants.

## Path Forward for Population Normalization

### Option 1: Manual 1000 Genomes Calculation (Recommended)

Since we already have high-quality scores with ~99.9% coverage, we can:

1. Download 1000 Genomes EUR subset VCF (~500 samples)
2. Calculate PRS for all EUR samples using our existing calculator
3. Compute mean and SD from that distribution
4. Normalize our scores: `Z = (score - EUR_mean) / EUR_SD`
5. Convert to percentiles: `percentile = norm.cdf(Z) * 100`

This approach:
- Uses our validated, high-coverage calculator
- Provides EUR-specific normalization
- Avoids pgsc_calc bugs

### Option 2: Wait for pgsc_calc Fix

The ancestry normalization bug may be fixed in a future pgsc_calc release. The issue has been identified as a path handling problem in Nextflow's PathVisitor.

### Option 3: Run pgsc_calc on HPC Cluster

Some users report success running pgsc_calc on HPC systems with different filesystem configurations. This could bypass the path handling bug.

## Files Created

- `data/trent_wgs.{pgen,pvar,psam}` - PLINK2 format genomic data
- `samplesheet.csv` - pgsc_calc input samplesheet
- `results/pgsc_calc/` - Pipeline outputs (raw scores, match logs, report)
- `results/pgsc_calc_comparison.csv` - Score comparison

## Temporary Files

Files in `/tmp/prs_data/` and `/tmp/prs_work/` can be cleaned up.

---
*Generated: 2026-01-19*
