#!/usr/bin/env python3
"""
Polygenic Risk Score Calculator v2.0
Handles variant-only VCF files with reference genome lookup for missing positions

Improvements over v18:
- Support for both GRCh37 and GRCh38 genome builds
- Command-line argument parsing
- Better error handling and validation
- Population/ancestry warnings
- Type hints throughout

Usage:
    python pgs_calculator_v2.py \
        --pgs-id PGS004237 \
        --vcf your_file.vcf \
        --reference /path/to/reference.fasta \
        --build GRCh38
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import requests


@dataclass
class PGSVariant:
    """Represents a variant from the PGS scoring file."""

    effect_allele: str
    other_allele: str | None
    weight: float
    rsid: str | None


@dataclass
class VCFVariant:
    """Represents a variant from the VCF file."""

    ref: str
    alt: str
    genotype: list[int]
    rsid: str | None


@dataclass
class ScoredVariant:
    """Result of scoring a single variant."""

    chrom: str
    pos: int
    rsid: str
    ref: str
    alt: str
    genotype: str
    effect_allele: str
    source: str  # 'vcf' or 'ref_lookup'
    effect_is_ref: bool
    dosage: int
    weight: float
    contribution: float


def download_pgs_file(pgs_id: str, output_file: str, build: str = "GRCh38") -> bool:
    """
    Download PGS Catalog scoring file.

    Args:
        pgs_id: PGS Catalog ID (e.g., 'PGS004237')
        output_file: Path to save the downloaded file
        build: Genome build ('GRCh37' or 'GRCh38')

    Returns:
        True if download successful, False otherwise
    """
    # Try harmonized version first
    url = (
        f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/"
        f"ScoringFiles/Harmonized/{pgs_id}_hmPOS_{build}.txt.gz"
    )

    print(f"Downloading {pgs_id} (harmonized {build}) from PGS Catalog...")
    print(f"URL: {url}")

    response = requests.get(url, timeout=60)

    if response.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(response.content)
        print(f"Downloaded to {output_file}")
        return True

    print(f"Error downloading harmonized version: HTTP {response.status_code}")
    print("Trying non-harmonized version...")

    # Try non-harmonized version
    url = (
        f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{pgs_id}/"
        f"ScoringFiles/{pgs_id}.txt.gz"
    )
    response = requests.get(url, timeout=60)

    if response.status_code == 200:
        with open(output_file, "wb") as f:
            f.write(response.content)
        print(f"Downloaded non-harmonized version to {output_file}")
        print(
            "WARNING: Non-harmonized file may not have positions in your target build!"
        )
        return True

    print(f"Error downloading: HTTP {response.status_code}")
    return False


def check_reference_genome(ref_file: str) -> str | None:
    """
    Check if reference genome file exists.

    Args:
        ref_file: Path to reference FASTA file

    Returns:
        Path to reference file if found, None otherwise
    """
    if os.path.exists(ref_file):
        print(f"Using reference genome: {ref_file}")
        # Check for index file
        if os.path.exists(ref_file + ".fai"):
            print(f"Index file found: {ref_file}.fai")
        return ref_file

    print(f"ERROR: Reference file not found: {ref_file}")
    return None


def load_reference_sequences(
    ref_file: str, chromosomes_needed: set[str]
) -> dict[str, str]:
    """
    Load reference sequences for specified chromosomes.

    Args:
        ref_file: Path to reference FASTA file
        chromosomes_needed: Set of chromosome names to load

    Returns:
        Dictionary mapping chromosome names to sequences

    Note:
        This loads entire chromosomes into memory - requires ~3GB RAM
    """
    sequences: dict[str, str] = {}

    print("Loading reference sequences...")
    print(f"Chromosomes needed: {sorted(chromosomes_needed)}")

    opener: Any = gzip.open if ref_file.endswith(".gz") else open

    current_chr: str | None = None
    current_seq: list[str] = []

    with opener(ref_file, "rt") as f:
        for line in f:
            if line.startswith(">"):
                # Save previous chromosome
                if current_chr and current_chr in chromosomes_needed:
                    sequences[current_chr] = "".join(current_seq)
                    print(f"  Loaded chr{current_chr}: {len(sequences[current_chr]):,} bp")

                # Start new chromosome
                # Handle different formats: >1, >chr1, >chrM, etc.
                chr_line = line.strip().split()[0][1:]  # Remove '>' and take first field
                if chr_line.startswith("chr"):
                    current_chr = chr_line[3:]
                else:
                    current_chr = chr_line

                current_seq = []

                # Skip if we don't need this chromosome
                if current_chr not in chromosomes_needed:
                    current_chr = None

            elif current_chr:
                current_seq.append(line.strip().upper())

    # Save last chromosome
    if current_chr and current_chr in chromosomes_needed:
        sequences[current_chr] = "".join(current_seq)
        print(f"  Loaded chr{current_chr}: {len(sequences[current_chr]):,} bp")

    print(f"Loaded {len(sequences)} chromosome sequences")
    return sequences


def get_reference_allele(sequences: dict[str, str], chrom: str, pos: int) -> str | None:
    """
    Get reference allele at specific position.

    Args:
        sequences: Dictionary of chromosome sequences
        chrom: Chromosome name
        pos: 1-based position (standard VCF/genomics convention)

    Returns:
        Reference allele or None if not found
    """
    if chrom not in sequences:
        return None

    # Convert to 0-based for Python string indexing
    idx = pos - 1

    if idx < 0 or idx >= len(sequences[chrom]):
        return None

    return sequences[chrom][idx]


def parse_pgs_file(pgs_file: str) -> dict[tuple[str, int], PGSVariant]:
    """
    Parse PGS scoring file and extract variant information.

    Args:
        pgs_file: Path to PGS scoring file

    Returns:
        Dictionary mapping (chr, pos) to PGSVariant objects
    """
    variants: dict[tuple[str, int], PGSVariant] = {}
    skipped_no_pos = 0

    opener: Any = gzip.open if pgs_file.endswith(".gz") else open

    with opener(pgs_file, "rt") as f:
        header: list[str] | None = None
        for line in f:
            if line.startswith("#"):
                continue

            if header is None:
                header = line.strip().split("\t")
                # Find column indices - use harmonized columns if available
                chr_idx = None
                pos_idx = None

                # Prefer harmonized columns (hm_chr, hm_pos)
                if "hm_chr" in header:
                    chr_idx = header.index("hm_chr")
                elif "chr_name" in header:
                    chr_idx = header.index("chr_name")

                if "hm_pos" in header:
                    pos_idx = header.index("hm_pos")
                elif "chr_position" in header:
                    pos_idx = header.index("chr_position")

                effect_idx = header.index("effect_allele")
                other_idx = header.index("other_allele") if "other_allele" in header else None
                weight_idx = header.index("effect_weight")

                # Try harmonized rsID first
                rsid_idx = None
                if "hm_rsID" in header:
                    rsid_idx = header.index("hm_rsID")
                elif "rsID" in header:
                    rsid_idx = header.index("rsID")

                continue

            # Use rstrip('\n') instead of strip() to preserve trailing tabs for empty columns
            fields = line.rstrip('\n').split("\t")

            # Skip if no position information or insufficient fields
            if chr_idx is None or pos_idx is None:
                skipped_no_pos += 1
                continue
            if len(fields) <= pos_idx or not fields[pos_idx] or fields[pos_idx] == "":
                skipped_no_pos += 1
                continue

            chrom = fields[chr_idx].replace("chr", "")
            pos = int(fields[pos_idx])
            effect_allele = fields[effect_idx]
            other_allele = fields[other_idx] if other_idx and len(fields) > other_idx and fields[other_idx] else None
            weight = float(fields[weight_idx])
            rsid = fields[rsid_idx] if rsid_idx and len(fields) > rsid_idx and fields[rsid_idx] else None

            variants[(chrom, pos)] = PGSVariant(
                effect_allele=effect_allele,
                other_allele=other_allele,
                weight=weight,
                rsid=rsid,
            )

    print(f"Parsed {len(variants)} variants from PGS file")
    if skipped_no_pos > 0:
        print(f"Skipped {skipped_no_pos} variants without position information")

    return variants


def load_vcf_variants(vcf_file: str) -> dict[tuple[str, int], VCFVariant]:
    """
    Load ALL variants from VCF into memory.

    Args:
        vcf_file: Path to VCF file

    Returns:
        Dictionary mapping (chr, pos) to VCFVariant objects
    """
    vcf_variants: dict[tuple[str, int], VCFVariant] = {}
    skipped = 0

    opener: Any = gzip.open if vcf_file.endswith(".gz") else open

    print("Loading all variants from VCF...")

    with opener(vcf_file, "rt") as f:
        line_count = 0
        for line in f:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")

            # Skip malformed lines
            if len(fields) < 10:
                continue

            line_count += 1
            if line_count % 1_000_000 == 0:
                print(f"  Loaded {line_count:,} variants...")

            chrom = fields[0].replace("chr", "")
            pos = int(fields[1])
            vcf_rsid = fields[2] if fields[2] != "." else None
            ref = fields[3]
            alt = fields[4]

            # Get genotype (assuming single sample, last column)
            gt_info = fields[9].split(":")[0]

            # Parse genotype
            if "/" in gt_info:
                alleles = gt_info.split("/")
            elif "|" in gt_info:
                alleles = gt_info.split("|")
            else:
                skipped += 1
                continue

            # Skip if missing genotype
            if "." in alleles:
                skipped += 1
                continue

            vcf_variants[(chrom, pos)] = VCFVariant(
                ref=ref,
                alt=alt,
                genotype=[int(a) for a in alleles],
                rsid=vcf_rsid,
            )

    print(f"Loaded {len(vcf_variants):,} variants from VCF")
    if skipped > 0:
        print(f"Skipped {skipped} variants with missing genotypes")

    return vcf_variants


def calculate_pgs(
    pgs_variants: dict[tuple[str, int], PGSVariant],
    vcf_variants: dict[tuple[str, int], VCFVariant],
    ref_sequences: dict[str, str],
) -> tuple[float, list[ScoredVariant]]:
    """
    Calculate PGS with proper handling of variant-only VCF and reference lookup.

    Args:
        pgs_variants: Dictionary of PGS variants
        vcf_variants: Dictionary of VCF variants
        ref_sequences: Dictionary of reference sequences

    Returns:
        Tuple of (total_score, list of ScoredVariant results)
    """
    total_score = 0.0
    scored_variants = 0

    # Track different scenarios
    found_in_vcf = 0
    looked_up_in_ref = 0
    ref_lookup_failed = 0
    allele_mismatch = 0

    ref_effect_in_vcf = 0
    alt_effect_in_vcf = 0
    ref_effect_homref = 0
    alt_effect_homref = 0

    results: list[ScoredVariant] = []

    for (chrom, pos), pgs_info in pgs_variants.items():
        effect_allele = pgs_info.effect_allele
        other_allele = pgs_info.other_allele
        weight = pgs_info.weight

        # Check if this position is in the VCF
        if (chrom, pos) in vcf_variants:
            # Position found in VCF - we have actual genotype data
            found_in_vcf += 1

            vcf_info = vcf_variants[(chrom, pos)]
            ref = vcf_info.ref
            alt = vcf_info.alt
            gt = vcf_info.genotype
            vcf_rsid = vcf_info.rsid

            # Use VCF rsID if available, otherwise use PGS rsID
            rsid_to_use = vcf_rsid if vcf_rsid else (pgs_info.rsid or ".")

            # Determine which allele is the effect allele
            effect_is_ref = effect_allele == ref
            effect_is_alt = effect_allele == alt

            # Check for allele match
            if not effect_is_ref and not effect_is_alt:
                allele_mismatch += 1
                continue

            # Calculate dosage of effect allele
            if effect_is_alt:
                dosage = sum(gt)  # Count ALT alleles
                alt_effect_in_vcf += 1
            else:
                dosage = 2 - sum(gt)  # Count REF alleles
                ref_effect_in_vcf += 1

            # Determine genotype string
            if gt == [0, 0]:
                genotype = f"{ref}/{ref}"
            elif gt in ([0, 1], [1, 0]):
                genotype = f"{ref}/{alt}"
            else:
                genotype = f"{alt}/{alt}"

            contribution = dosage * weight
            total_score += contribution
            scored_variants += 1

            results.append(
                ScoredVariant(
                    chrom=chrom,
                    pos=pos,
                    rsid=rsid_to_use,
                    ref=ref,
                    alt=alt,
                    genotype=genotype,
                    effect_allele=effect_allele,
                    source="vcf",
                    effect_is_ref=effect_is_ref,
                    dosage=dosage,
                    weight=weight,
                    contribution=contribution,
                )
            )

        else:
            # Position NOT in VCF - look up reference allele
            looked_up_in_ref += 1

            ref_allele = get_reference_allele(ref_sequences, chrom, pos)

            if ref_allele is None:
                ref_lookup_failed += 1
                continue

            # Determine if effect allele matches reference
            effect_is_ref = effect_allele == ref_allele

            if effect_is_ref:
                # Effect allele IS the reference, and we're homozygous reference
                dosage = 2
                ref_effect_homref += 1
            else:
                # Effect allele is NOT the reference (must be ALT)
                # We're homozygous reference, so we have ZERO copies of effect allele
                dosage = 0
                alt_effect_homref += 1

            genotype = f"{ref_allele}/{ref_allele}"

            contribution = dosage * weight
            total_score += contribution
            scored_variants += 1

            results.append(
                ScoredVariant(
                    chrom=chrom,
                    pos=pos,
                    rsid=pgs_info.rsid or ".",
                    ref=ref_allele,
                    alt=other_allele or "?",
                    genotype=genotype,
                    effect_allele=effect_allele,
                    source="ref_lookup",
                    effect_is_ref=effect_is_ref,
                    dosage=dosage,
                    weight=weight,
                    contribution=contribution,
                )
            )

    print("\n=== PGS Calculation Summary ===")
    print(f"Total variants in PGS file: {len(pgs_variants):,}")
    print(f"\nVariant sources:")
    print(f"  Found in VCF: {found_in_vcf:,}")
    print(f"  Looked up in reference: {looked_up_in_ref:,}")
    print(f"  Reference lookup failed: {ref_lookup_failed:,}")
    print(f"  Allele mismatches (excluded): {allele_mismatch:,}")
    print(f"\nEffect allele breakdown:")
    print(f"  From VCF:")
    print(f"    Effect = ALT: {alt_effect_in_vcf:,}")
    print(f"    Effect = REF: {ref_effect_in_vcf:,}")
    print(f"  From reference lookup (homozygous reference):")
    print(f"    Effect = REF (dosage=2): {ref_effect_homref:,}")
    print(f"    Effect = ALT (dosage=0): {alt_effect_homref:,}")
    print(f"\nSuccessfully scored: {scored_variants:,}")
    print(f"Coverage: {100 * scored_variants / len(pgs_variants):.1f}%")
    print(f"\n{'=' * 50}")
    print(f"YOUR POLYGENIC SCORE: {total_score:.4f}")
    print(f"{'=' * 50}")

    return total_score, results


def save_results(results: list[ScoredVariant], output_file: str) -> None:
    """Save detailed results to a tab-separated file."""
    with open(output_file, "w") as f:
        f.write(
            "chr\tpos\trsID\tref\talt\tgenotype\teffect_allele\t"
            "source\teffect_is_ref\tdosage\tweight\tcontribution\n"
        )
        for r in results:
            f.write(
                f"{r.chrom}\t{r.pos}\t{r.rsid}\t{r.ref}\t{r.alt}\t"
                f"{r.genotype}\t{r.effect_allele}\t{r.source}\t{r.effect_is_ref}\t"
                f"{r.dosage}\t{r.weight:.6f}\t{r.contribution:.6f}\n"
            )
    print(f"\nDetailed results saved to {output_file}")


def print_ancestry_warning() -> None:
    """Print warning about ancestry and PGS interpretation."""
    print("\n" + "=" * 60)
    print("IMPORTANT: PGS Ancestry Considerations")
    print("=" * 60)
    print(
        """
Most PGS scores have been developed using European ancestry populations.
If your genetic ancestry differs, the score may be:
  - Less predictive of actual risk
  - Not directly comparable to published reference distributions
  
To understand your score in context, consider:
  1. The ancestry composition of the discovery GWAS
  2. The ancestry of validation cohorts
  3. Your own genetic ancestry

For ancestry-adjusted scores, consider using the PGS Catalog Calculator
with --run_ancestry flag, which uses reference populations (HGDP+1kGP).
"""
    )


def parse_args() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate Polygenic Risk Score from VCF file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with GRCh38
  python pgs_calculator_v2.py --pgs-id PGS004237 --vcf my.vcf --reference hg38.fa

  # Specify genome build explicitly
  python pgs_calculator_v2.py --pgs-id PGS004237 --vcf my.vcf --reference hg38.fa --build GRCh38

  # Use pre-downloaded PGS file
  python pgs_calculator_v2.py --pgs-file PGS004237_hmPOS_GRCh38.txt.gz --vcf my.vcf --reference hg38.fa
        """,
    )

    parser.add_argument(
        "--pgs-id",
        help="PGS Catalog ID (e.g., PGS004237). Will download if --pgs-file not provided.",
    )
    parser.add_argument(
        "--pgs-file",
        help="Path to pre-downloaded PGS scoring file (skips download)",
    )
    parser.add_argument(
        "--vcf",
        required=True,
        help="Path to VCF file containing your variants",
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference genome FASTA file",
    )
    parser.add_argument(
        "--build",
        choices=["GRCh37", "GRCh38"],
        default="GRCh38",
        help="Genome build (default: GRCh38)",
    )
    parser.add_argument(
        "--output",
        help="Output file prefix (default: <PGS_ID>_results)",
    )
    parser.add_argument(
        "--no-ancestry-warning",
        action="store_true",
        help="Suppress ancestry warning message",
    )

    args = parser.parse_args()

    # Validate that either pgs-id or pgs-file is provided
    if not args.pgs_id and not args.pgs_file:
        parser.error("Either --pgs-id or --pgs-file must be provided")

    return args


def main() -> None:
    """Main entry point."""
    args = parse_args()

    # Determine PGS file path
    if args.pgs_file:
        pgs_file = args.pgs_file
        pgs_id = Path(args.pgs_file).stem.split("_")[0]
    else:
        pgs_id = args.pgs_id
        pgs_file = f"{pgs_id}_hmPOS_{args.build}.txt.gz"

    # Step 1: Get PGS file
    print("=" * 60)
    print("STEP 1: Get PGS scoring file")
    print("=" * 60)

    if args.pgs_file and os.path.exists(args.pgs_file):
        print(f"Using provided PGS file: {args.pgs_file}")
    else:
        if not download_pgs_file(pgs_id, pgs_file, args.build):
            print("Failed to download PGS file")
            sys.exit(1)

    # Step 2: Parse PGS file
    print("\n" + "=" * 60)
    print("STEP 2: Parse PGS scoring file")
    print("=" * 60)
    pgs_variants = parse_pgs_file(pgs_file)

    # Get list of chromosomes we need
    chromosomes_needed = {chrom for chrom, pos in pgs_variants.keys()}

    # Step 3: Check for reference genome
    print("\n" + "=" * 60)
    print("STEP 3: Load reference genome")
    print("=" * 60)
    ref_file = check_reference_genome(args.reference)
    if not ref_file:
        print("Cannot proceed without reference genome")
        sys.exit(1)

    ref_sequences = load_reference_sequences(ref_file, chromosomes_needed)

    # Step 4: Load ALL variants from VCF
    print("\n" + "=" * 60)
    print("STEP 4: Load VCF variants")
    print("=" * 60)
    vcf_variants = load_vcf_variants(args.vcf)

    # Step 5: Calculate PGS
    print("\n" + "=" * 60)
    print("STEP 5: Calculate polygenic score")
    print("=" * 60)
    score, results = calculate_pgs(pgs_variants, vcf_variants, ref_sequences)

    # Step 6: Save detailed results
    output_prefix = args.output or f"{pgs_id}_results"
    output_file = f"{output_prefix}.txt"
    save_results(results, output_file)

    # Print ancestry warning
    if not args.no_ancestry_warning:
        print_ancestry_warning()


if __name__ == "__main__":
    main()
