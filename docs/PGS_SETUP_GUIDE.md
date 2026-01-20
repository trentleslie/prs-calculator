# PGS Calculator Setup Guide

## Overview

This guide covers three goals from the Slack conversation:

1. **Calculate your own PRS scores** with your VCF file
2. **Build a Replit app** for users to upload VCF files and calculate PRS
3. **Determine demographic/ancestry information** from your 23andMe-formatted file

---

## Goal 1: Calculate Your PRS Scores

### Prerequisites

1. **Your VCF file** - Likely from Arivale/Nucleus sequencing
2. **Reference genome FASTA** - Must match your VCF's build:
   - GRCh37/hg19: `homo_sapiens_assembly19.fasta` (Lee's setup)
   - GRCh38/hg38: Check the Phenome core (Gwen posted link in DS channel)

### Determining Your VCF's Genome Build

Before running, verify which build your VCF uses:

```bash
# Check VCF header for build info
grep -m 20 "^#" your_file.vcf | grep -i "reference\|assembly\|build"

# Or compare a known position - chr1 length differs between builds:
# GRCh37: 249,250,621 bp
# GRCh38: 248,956,422 bp
```

Arivale VCFs were typically GRCh37/hg19. Nucleus sequencing may be GRCh38.

### Running the Calculator

```bash
# For GRCh38 (if your VCF is aligned to hg38)
python pgs_calculator_v2.py \
    --pgs-id PGS004237 \
    --vcf /path/to/your.vcf \
    --reference /path/to/GRCh38_reference.fasta \
    --build GRCh38

# For GRCh37 (if your VCF is aligned to hg19)
python pgs_calculator_v2.py \
    --pgs-id PGS004237 \
    --vcf /path/to/your.vcf \
    --reference /path/to/homo_sapiens_assembly19.fasta \
    --build GRCh37
```

### Finding PGS Scores to Explore

Visit [PGS Catalog](https://www.pgscatalog.org/) to find scores for traits you're interested in.

**Important**: As Lee mentioned, choose scores that handle Linkage Disequilibrium (LD):

| LD-Aware Methods (GOOD) | Methods Without LD (AVOID) |
|------------------------|---------------------------|
| LDpred, LDpred2 | Simple sum of GWAS betas |
| PRS-CS | Unweighted counts |
| lassosum | P+T without clumping |
| SBayesR | |
| PRSice-2 (clumping) | |

Check the PGS Catalog entry's "Method" field to see which approach was used.

### Interpreting Your Score

- Raw PGS scores are typically centered around 0
- Positive scores generally indicate higher genetic predisposition
- **Ancestry matters**: Most PGS were developed from European populations
- Compare to population distributions when available

---

## Goal 2: Replit App Development

### App Concept

A web application where users can:
1. Upload their VCF file (or drag-and-drop)
2. Select a trait/PGS from a curated list
3. Get their calculated PGS with interpretation

### Replit Prompt

Here's a prompt to get started:

```
Build a Flask/Python web application for calculating Polygenic Risk Scores (PGS).

Core Features:
1. File upload for VCF files (support .vcf and .vcf.gz)
2. Dropdown to select from popular PGS Catalog entries
3. Backend that:
   - Downloads PGS scoring files from PGS Catalog FTP
   - Parses the VCF to extract genotypes
   - Calculates the polygenic score
   - Returns results with basic interpretation

Technical Requirements:
- Use Flask for the web framework
- Store uploaded VCFs temporarily (delete after processing)
- Cache downloaded PGS files to avoid repeated downloads
- Show progress indicator during calculation
- Include ancestry warning/disclaimer

UI/UX:
- Simple, clean interface
- Clear explanation of what PGS is
- Disclaimer about limitations (ancestry, not medical advice)
- Download results as CSV

Security Considerations:
- Validate uploaded files
- Size limits on uploads
- No persistent storage of genetic data
- HTTPS only

Reference Implementation:
- Use the algorithm from proper_pgs_calculator_v18.py
- PGS Catalog FTP: https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
- Download harmonized files: {PGS_ID}_hmPOS_GRCh38.txt.gz

Note: For a Replit app, we'll use a pre-bundled minimal reference 
or rsID-based matching to avoid large reference genome files.
```

### Architecture Considerations

1. **Reference Genome Challenge**: Full reference genomes are ~3GB
   - Option A: Use rsID matching instead of position-based (requires VCF with rsIDs)
   - Option B: Pre-compute reference alleles for common PGS variants
   - Option C: Use a cloud storage bucket for reference data

2. **Performance**: 
   - Cache PGS scoring files
   - Consider WebAssembly for client-side processing
   - Use streaming for large VCF files

---

## Goal 3: Ancestry/Demographics from 23andMe File

### What the 23andMe Format Contains

The "23AndMe Formatted File" from Nucleus is a simplified extraction containing:
- rsID (variant identifier)
- Chromosome
- Position
- Genotype (e.g., AA, AG, GG)

### Ancestry Analysis Options

#### Option A: Third-Party Services
Upload to services that accept 23andMe format:
- **GEDmatch** - Free ancestry analysis
- **Promethease** - Health/trait analysis (~$12)
- **DNA.Land** - Free ancestry composition

#### Option B: DIY Ancestry with PCA

You can estimate ancestry using Principal Component Analysis with reference populations:

```python
"""
Ancestry estimation using PCA with 1000 Genomes reference panel

This approach:
1. Extracts AIMs (Ancestry Informative Markers) from your file
2. Projects onto PCA space computed from 1000 Genomes populations
3. Estimates admixture proportions
"""

# Required: 
# - 1000 Genomes VCF (or pre-computed PCA loadings)
# - PLINK for PCA
# - Your 23andMe file converted to PLINK format

# The PGS Catalog Calculator supports ancestry analysis:
# https://pgsc-calc.readthedocs.io/en/latest/
# Use: --run_ancestry flag with reference panel
```

#### Option C: Simpler Marker-Based Approach

Check specific ancestry-informative SNPs:

```python
"""
Check specific markers associated with population ancestry.
This is a simplified approach - not as accurate as full PCA.
"""

# Example AIMs (Ancestry Informative Markers)
ANCESTRY_MARKERS = {
    # African ancestry markers
    'rs2814778': {'CC': 'African', 'TT': 'non-African'},  # Duffy null
    
    # European ancestry markers  
    'rs16891982': {'GG': 'European', 'CC': 'non-European'},  # SLC45A2
    
    # East Asian markers
    'rs3827760': {'AA': 'East Asian', 'GG': 'non-Asian'},  # EDAR
    
    # Native American markers
    # ...
}

def check_ancestry_markers(genotypes: dict) -> dict:
    """
    Check genotypes at ancestry-informative markers.
    Returns dict of marker results.
    """
    results = {}
    for rs, associations in ANCESTRY_MARKERS.items():
        if rs in genotypes:
            gt = genotypes[rs]
            results[rs] = associations.get(gt, 'Unknown')
    return results
```

### Recommended Approach for Ancestry

Given you have full WGS data, the best approach would be:

1. **Use the PGS Catalog Calculator** with `--run_ancestry` flag
   - Downloads HGDP+1kGP reference panel
   - Performs proper PCA-based ancestry estimation
   - Provides population similarity scores

2. **Check Nucleus directly** - They may have ancestry analysis available

3. **Use GEDmatch** with your 23andMe-formatted file for quick results

---

## Quick Reference: File Locations

Based on the Slack conversation:

| Resource | Location |
|----------|----------|
| GRCh38 Reference | Phenome core (link in DS channel) |
| GRCh37 Reference | `homo_sapiens_assembly19.fasta` |
| Lee's VCF | `LR_full_variant_file.vcf` |
| PGS Catalog | https://www.pgscatalog.org |
| Score FTP | https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/ |

---

## Troubleshooting

### VCF Format Issues

Your VCF should have 10 columns (standard format). If it differs:
- Check if it's a multi-sample VCF
- Verify column delimiters (should be tabs)
- Check for extra header lines

### Score Interpretation

- **Coverage < 90%**: May indicate build mismatch or missing variants
- **Very extreme scores**: Double-check allele orientation
- **Score = 0**: Usually indicates a problem with variant matching

### Memory Issues

Loading the reference genome requires ~3GB RAM. If memory-constrained:
- Use a machine with more RAM
- Consider chunked processing
- Use rsID-based matching instead
