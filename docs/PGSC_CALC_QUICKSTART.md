# pgsc_calc Quickstart Guide

This guide shows you how to run pgsc_calc on a **single VCF file** from whole genome sequencing. This is the most common use case that isn't well documented.

---

## Prerequisites

### Required Software

| Tool | Purpose | Installation |
|------|---------|--------------|
| **Nextflow** | Workflow engine | `curl -s https://get.nextflow.io \| bash` |
| **Conda** or **Docker** | Environment management | Conda recommended for local runs |
| **PLINK2** | VCF conversion | `conda install -c bioconda plink2` |
| **bcftools** | VCF utilities | `conda install -c bioconda bcftools` |

### Verify Installation

```bash
nextflow -version    # Should be 23.04+
plink2 --version     # Should be 2.0+
bcftools --version   # Should be 1.15+
```

---

## Quick Answer: Running on a Single VCF

### Step 1: Determine Your VCF's Reference Build

This is **critical** - mismatched builds will produce garbage results.

```bash
# Check VCF header for reference info
bcftools view -h your_file.vcf.gz | grep -i "reference"

# Alternative: check contig lengths
# chr1 length: 248,956,422 bp (GRCh38) vs 249,250,621 bp (GRCh37)
bcftools view -h your_file.vcf.gz | grep "##contig=<ID=chr1"
```

**Our Nucleus data is GRCh38** (Illumina DRAGEN pipeline).

---

### Step 2: Convert VCF to PLINK2 Format

pgsc_calc requires PLINK2 pgen/pvar/psam format, not VCF.

```bash
# Create output directory
mkdir -p pgsc_calc_input

# Convert VCF to PLINK2 format
plink2 --vcf your_file.vcf.gz \
    --make-pgen \
    --out pgsc_calc_input/your_sample \
    --allow-extra-chr \
    --vcf-half-call m \
    --threads 4

# Verify output (should have 3 files)
ls pgsc_calc_input/
# your_sample.pgen  your_sample.pvar  your_sample.psam
```

**Note**: The `--vcf-half-call m` flag treats half-calls (e.g., `0/.`) as missing, which is the safest default.

---

### Step 3: Create Samplesheet

pgsc_calc requires a CSV samplesheet describing your input files.

```bash
# Create samplesheet.csv
cat > samplesheet.csv << 'EOF'
sampleset,path_prefix,chrom,format
trent,/absolute/path/to/pgsc_calc_input/your_sample,,pfile
EOF
```

**Important:**
- Use **absolute paths** (pgsc_calc runs in a container/conda env)
- Leave `chrom` empty (we have whole genome, not split by chromosome)
- Format is `pfile` for PLINK2 pgen/pvar/psam format

---

### Step 4: Run pgsc_calc

```bash
# Basic run with specific PGS scores
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input samplesheet.csv \
    --pgs_id PGS002308,PGS004034 \
    --target_build GRCh38 \
    -r v2.0.0-alpha

# With ancestry normalization (recommended)
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input samplesheet.csv \
    --pgs_id PGS002308,PGS004034 \
    --target_build GRCh38 \
    --run_ancestry ALL \
    --ancestry_runtype direct \
    -r v2.0.0-alpha
```

**Parameter explanation:**

| Parameter | Purpose |
|-----------|---------|
| `-profile conda` | Use conda for dependencies (or `docker`/`singularity`) |
| `--pgs_id` | Comma-separated PGS Catalog IDs to calculate |
| `--target_build` | Your VCF's reference build (GRCh37 or GRCh38) |
| `--run_ancestry ALL` | Enable ancestry normalization |
| `--ancestry_runtype direct` | Use built-in reference panel |
| `-r v2.0.0-alpha` | Pin to specific version (recommended) |

---

### Step 5: View Results

```bash
# Results are in results/ directory
ls results/

# Open interactive report
open results/trent/report.html  # macOS
xdg-open results/trent/report.html  # Linux

# View summary scores
zcat results/trent/trent_pgs.txt.gz | head -20

# View ancestry analysis
zcat results/trent/trent_popsimilarity.txt.gz
```

---

## GRCh37/hg19 Support

pgsc_calc **does** support GRCh37, but with important caveats.

### Using GRCh37 Data

```bash
# Same process, just specify GRCh37 build
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input samplesheet.csv \
    --pgs_id PGS002308 \
    --target_build GRCh37 \
    -r v2.0.0-alpha
```

### Important Limitations

| Issue | Details |
|-------|---------|
| **Score availability** | Some newer PGS scores are GRCh38-only |
| **No automatic liftover** | pgsc_calc expects harmonized scores for your build |
| **Reference panel** | HGDP+1kGP panel is GRCh38; ancestry normalization may be limited |

### Checking Score Availability

Not all PGS Catalog scores have GRCh37 harmonized versions:

```bash
# Check if GRCh37 version exists
curl -s "https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002308/ScoringFiles/Harmonized/" | grep GRCh37
```

If no GRCh37 version exists, you have two options:
1. Lift over your VCF to GRCh38 (recommended)
2. Use a custom script with manual coordinate conversion

---

## What About LD Correction?

This is a common point of confusion.

### Key Insight: pgsc_calc Does NOT Perform LD Correction

The LD correction happens **earlier in the pipeline**, not at scoring time:

```
1. GWAS identifies associations (raw, LD-confounded)
       ↓
2. PGS developers apply LD correction (LDpred, PRS-CS, etc.)
       ↓
3. PGS Catalog publishes the corrected weights
       ↓
4. pgsc_calc applies those pre-corrected weights to your genotypes
```

### Why This Matters

When you select an "LD-aware" PGS like:
- **PGS002308** (T2D, uses PRS-CSx)
- **PGS004034** (AD, uses LDpred2-auto)

The LD handling was done by the **score developers**, not pgsc_calc. The weights you download already account for LD structure.

### What pgsc_calc DOES Handle

| Step | What pgsc_calc does |
|------|---------------------|
| Variant matching | Maps your VCF variants to PGS scoring file |
| Allele alignment | Flips alleles if REF/ALT are swapped |
| Missing variants | Imputes from reference panel allele frequencies |
| Ancestry normalization | Adjusts scores relative to population distribution |

### What pgsc_calc Does NOT Do

- Apply LD correction to raw GWAS betas
- Clump variants
- Prune correlated variants

These are all handled upstream by the PGS score developers.

---

## Troubleshooting

### "No variants matched"

**Cause**: Build mismatch between VCF and scoring file

**Fix**: Verify your VCF build and use matching `--target_build`

```bash
# Check your VCF build
bcftools view -h your_file.vcf.gz | grep "##reference"
```

### "Can't find reference panel"

**Cause**: Ancestry normalization requires downloading reference data

**Fix**: Run without ancestry normalization first, or pre-download panel

```bash
# Run without ancestry normalization
nextflow run pgscatalog/pgsc_calc \
    --input samplesheet.csv \
    --pgs_id PGS002308 \
    --target_build GRCh38 \
    -profile conda
```

### Out of Memory

**Cause**: Large PGS scores (1M+ variants) need significant RAM

**Fix**: Increase memory allocation

```bash
# For Nextflow with increased memory
nextflow run pgscatalog/pgsc_calc \
    ... \
    -process.memory 16.GB
```

### Conda Environment Issues

**Cause**: Nextflow can't create conda environments

**Fix**: Use Docker instead, or pre-create environments

```bash
# Use Docker instead of conda
nextflow run pgscatalog/pgsc_calc \
    ... \
    -profile docker
```

---

## Example: Full Workflow

Here's a complete example using our Nucleus Genomics data:

```bash
# 1. Set up directory
cd Analysis/PRS
mkdir -p pgsc_calc_input results

# 2. Convert VCF to PLINK2
plink2 --vcf ../../Raw_Data/VCF/*.vcf.gz \
    --make-pgen \
    --out pgsc_calc_input/trent \
    --allow-extra-chr \
    --vcf-half-call m

# 3. Create samplesheet
cat > samplesheet.csv << EOF
sampleset,path_prefix,chrom,format
trent,$(pwd)/pgsc_calc_input/trent,,pfile
EOF

# 4. Run pgsc_calc with multiple scores
nextflow run pgscatalog/pgsc_calc \
    -profile conda \
    --input samplesheet.csv \
    --pgs_id PGS002308,PGS004034,PGS000027,PGS004237 \
    --target_build GRCh38 \
    --run_ancestry ALL \
    --ancestry_runtype direct \
    --outdir results/pgsc_calc \
    -r v2.0.0-alpha

# 5. View results
open results/pgsc_calc/trent/report.html
```

---

## Custom Script Alternative

For learning, debugging, or GRCh37 data, our custom Python script is available:

```bash
# Quick single-score calculation
uv run python src/pgs_calculator.py \
    --vcf ../../Raw_Data/VCF/*.vcf.gz \
    --pgs-file scores/PGS004237_score.txt \
    --output results/
```

**Advantages of custom script:**
- Runs in seconds (vs minutes for pgsc_calc)
- Shows every variant's contribution
- Works offline
- More transparent for debugging

**Disadvantages:**
- No ancestry normalization
- May have build compatibility issues (designed for GRCh37)
- Less robust handling of edge cases

See [PRS_METHODOLOGY_COMPARISON.md](PRS_METHODOLOGY_COMPARISON.md) for detailed comparison.

---

## Resources

- [pgsc_calc Documentation](https://pgsc-calc.readthedocs.io/)
- [PGS Catalog](https://www.pgscatalog.org/)
- [PLINK2 Documentation](https://www.cog-genomics.org/plink/2.0/)
- [PRS_METHODOLOGY_COMPARISON.md](PRS_METHODOLOGY_COMPARISON.md) - Our methodology comparison
