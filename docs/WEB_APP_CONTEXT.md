# Web App Integration Context

This document provides context for an LLM agent to build a web application that incorporates PRS (Polygenic Risk Score) calculation and visualization.

## Overview

This repository contains a complete PRS analysis workflow:
1. Download scoring files from PGS Catalog
2. Calculate raw PRS from VCF genotype data
3. Normalize against reference population (ancestry adjustment)
4. Visualize and interpret results

## Data Schemas

### Input: VCF File

Standard VCF format (v4.2+) with variants called against GRCh38/hg38.

```
##fileformat=VCFv4.2
##reference=GRCh38
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE
chr1    12345   rs123   A       G       99      PASS    .       GT      0/1
```

**Required columns:** CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, [SAMPLE]
**Required FORMAT field:** GT (genotype: 0/0, 0/1, 1/0, or 1/1)

### Input: PGS Scoring File

TSV format from [PGS Catalog](https://www.pgscatalog.org/).

| Column | Description | Example |
|--------|-------------|---------|
| `chr_name` | Chromosome (1-22, X, Y) | `1` |
| `chr_position` | Position in GRCh38 | `12345` |
| `effect_allele` | Allele that increases trait | `G` |
| `other_allele` | Reference/non-effect allele | `A` |
| `effect_weight` | Beta or log-odds coefficient | `0.0234` |

### Output: PRS Results JSON

```json
{
  "sample_id": "trent",
  "analysis_date": "2026-01-19",
  "reference_build": "GRCh38",
  "scores": [
    {
      "pgs_id": "PGS002308",
      "trait": "Type 2 Diabetes",
      "raw_score": 0.43638,
      "z_score": 2.76,
      "percentile": 99.7,
      "reference_population": "EUR",
      "variants_matched": 1045632,
      "variants_total": 1259754,
      "match_rate": 0.83,
      "interpretation": "Very High Genetic Risk"
    },
    {
      "pgs_id": "PGS004034",
      "trait": "Alzheimer's Disease",
      "raw_score": 1.67769,
      "z_score": 0.49,
      "percentile": 64.0,
      "reference_population": "EUR",
      "interpretation": "Average Genetic Risk"
    },
    {
      "pgs_id": "PGS000027",
      "trait": "Body Mass Index",
      "raw_score": 17.70405,
      "z_score": 0.06,
      "percentile": 51.6,
      "reference_population": "EUR",
      "interpretation": "Average Genetic Risk"
    },
    {
      "pgs_id": "PGS004237",
      "trait": "Coronary Artery Disease",
      "raw_score": -0.25798,
      "z_score": -1.47,
      "percentile": 7.2,
      "reference_population": "EUR",
      "interpretation": "Low Genetic Risk"
    }
  ],
  "ancestry": {
    "method": "PCA + Random Forest",
    "most_similar": "EUR",
    "confidence": "high",
    "note": "100% probability assigned to EUR"
  }
}
```

## Key Algorithms

### 1. PRS Calculation

```python
def calculate_prs(vcf_genotypes: dict, pgs_weights: list) -> float:
    """
    Calculate raw Polygenic Risk Score.

    Args:
        vcf_genotypes: Dict mapping chr:pos to genotype (0/0, 0/1, 1/1)
        pgs_weights: List of dicts with chr_position, effect_allele, effect_weight

    Returns:
        Raw PRS score (sum of dosage * weight)
    """
    score = 0.0
    matched = 0

    for variant in pgs_weights:
        key = f"{variant['chr_name']}:{variant['chr_position']}"
        if key in vcf_genotypes:
            genotype = vcf_genotypes[key]
            # Count copies of effect allele (dosage: 0, 1, or 2)
            dosage = count_effect_alleles(genotype, variant['effect_allele'], variant['other_allele'])
            score += dosage * variant['effect_weight']
            matched += 1

    return score, matched

def count_effect_alleles(genotype: str, effect: str, other: str) -> int:
    """
    Count effect allele copies from genotype string.

    Handles:
    - Standard: 0/0 (hom ref), 0/1 or 1/0 (het), 1/1 (hom alt)
    - Strand flips for ambiguous alleles (A/T, C/G)
    """
    alleles = genotype.replace('|', '/').split('/')
    count = 0
    for allele in alleles:
        if allele == '1':  # ALT allele
            count += 1
    return count
```

### 2. Ancestry Normalization (pgsc_calc method)

```python
def normalize_prs(raw_score: float, reference_panel: dict, ancestry: str) -> dict:
    """
    Normalize PRS against ancestry-matched reference population.

    Uses the pgsc_calc v2.0.0 methodology:
    1. Project sample onto reference panel PCA
    2. Classify ancestry using Random Forest
    3. Normalize within most similar population

    Args:
        raw_score: Raw PRS from calculate_prs()
        reference_panel: Dict with population -> array of scores
        ancestry: Most similar population from PCA (EUR, AFR, EAS, etc.)

    Returns:
        Dict with z_score and percentile
    """
    import scipy.stats as stats

    pop_scores = reference_panel[ancestry]
    pop_mean = np.mean(pop_scores)
    pop_std = np.std(pop_scores)

    z_score = (raw_score - pop_mean) / pop_std
    percentile = stats.norm.cdf(z_score) * 100

    return {
        'z_score': round(z_score, 2),
        'percentile': round(percentile, 1),
        'reference_population': ancestry
    }
```

### 3. Variant Matching Strategy

```python
def match_variants(vcf_variants: set, pgs_variants: list) -> list:
    """
    Match VCF variants to PGS scoring file variants.

    Matching hierarchy:
    1. Exact chr:pos:ref:alt match
    2. chr:pos match with allele flip (strand correction)
    3. chr:pos match with complementary alleles (A<>T, C<>G)

    Returns list of matched variants with dosage information.
    """
    COMPLEMENTS = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    matched = []
    for pgs_var in pgs_variants:
        pos_key = f"{pgs_var['chr_name']}:{pgs_var['chr_position']}"

        if pos_key in vcf_variants:
            vcf_var = vcf_variants[pos_key]

            # Check for allele match or flip
            if alleles_match(vcf_var, pgs_var, COMPLEMENTS):
                matched.append(create_match(vcf_var, pgs_var))

    return matched
```

## Web App Architecture Suggestions

### Backend (Python/FastAPI)

```python
# Recommended structure
app/
├── main.py              # FastAPI application
├── routers/
│   ├── upload.py        # VCF file upload
│   ├── calculate.py     # PRS calculation
│   └── results.py       # Results retrieval
├── services/
│   ├── vcf_parser.py    # cyvcf2-based parsing
│   ├── prs_calculator.py # Core calculation
│   ├── normalizer.py    # Ancestry normalization
│   └── pgs_catalog.py   # PGS file downloads
├── models/
│   ├── schemas.py       # Pydantic models
│   └── database.py      # SQLAlchemy models
└── workers/
    └── tasks.py         # Celery async tasks
```

**Key Dependencies:**
- `fastapi` - API framework
- `cyvcf2` - Fast VCF parsing (C-based)
- `pandas`, `numpy` - Data processing
- `scipy` - Statistical functions
- `celery` + `redis` - Async task queue (for large files)
- `sqlalchemy` - Database ORM

### Frontend Suggestions

- **Framework**: React, Next.js, or Vue
- **Visualization**: D3.js, Plotly.js, or Recharts
- **File Upload**: react-dropzone with progress indicator

### API Endpoints

```
POST /api/upload
  - Accepts: multipart/form-data with VCF file
  - Returns: { job_id: string, status: "processing" }

GET /api/jobs/{job_id}
  - Returns: { status: "pending"|"processing"|"complete"|"failed", progress: number }

GET /api/results/{job_id}
  - Returns: PRS Results JSON (see schema above)

GET /api/traits
  - Returns: List of available PGS scores with metadata

POST /api/calculate
  - Accepts: { job_id: string, pgs_ids: string[] }
  - Triggers calculation for specific scores
```

## Reference Implementation Files

| File | Purpose |
|------|---------|
| `src/pgs_calculator.py` | Core PRS calculation logic (variant matching, dosage counting) |
| `PRS_Analysis.ipynb` | Complete workflow demonstration with visualizations |
| `results/pgsc_calc/report.html` | Example interactive HTML report |
| `results/prs_normalized_ancestry.csv` | Output format example |
| `docs/PGS_SETUP_GUIDE.md` | Tool installation and setup |

## Clinical Interpretation Guide

| Percentile | Risk Level | Guidance |
|------------|------------|----------|
| >95% | Very High | Clinical genetics referral recommended |
| 75-95% | Elevated | Enhanced screening/prevention discussion |
| 25-75% | Average | Standard clinical care |
| 5-25% | Below Average | Lower genetic risk (lifestyle still matters) |
| <5% | Low | Possible genetic protection |

**Important Caveats:**
- PRS captures only genetic risk; environmental and lifestyle factors are equally important
- Ancestry matching is critical - scores derived from one population may not transfer well
- Clinical decisions should not be based solely on PRS
- New PGS scores are published regularly; consider version tracking

## Deployment Considerations

### VCF File Handling
- VCF files can be large (100MB-5GB for WGS)
- Consider streaming/chunked processing
- Validate file format before processing
- Support both `.vcf` and `.vcf.gz` formats

### PGS File Management
- Cache downloaded PGS files (they rarely change)
- Store version info for reproducibility
- Support multiple genome builds (hg19, hg38)

### Compute Resources
- Single-score calculation: ~30 seconds
- Full panel (4 scores): ~2 minutes
- Consider background workers for large batches

### Privacy/Security
- Genomic data is highly sensitive personal information
- Consider HIPAA compliance requirements
- Encrypt data at rest and in transit
- Implement proper access controls
- Consider processing locally vs. cloud upload
