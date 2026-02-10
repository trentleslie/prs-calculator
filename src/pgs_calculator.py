#!/usr/bin/env python3
"""
Proper PGS Calculator v33
Calculates BOTH adjusted (with reference lookup) and unadjusted (VCF only) scores
Outputs both to allow apples-to-apples comparison with reference panel

This version combines CLI argument support with dual-score calculation methodology.

Key features:
- Calculates ADJUSTED score: includes reference genome lookups (correct for WGS)
- Calculates UNADJUSTED score: VCF data only (matches pgsc_calc behavior)
- Supports --no-reference flag for VCF-only mode
- Command-line argument parsing for reusable execution
- Compatible with PRS_Methodology_Experiment.ipynb output format
"""

import argparse
import gzip
import os
import re
import sys
from pathlib import Path


def parse_vcf_header(vcf_file):
    """Parse VCF to get sample name"""
    opener = gzip.open if vcf_file.endswith('.gz') else open
    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                sample_name = headers[9] if len(headers) > 9 else 'SAMPLE'
                return sample_name
    return 'SAMPLE'


def load_pgs_catalog(pgs_file):
    """Load PGS Catalog file"""
    print(f"Loading PGS file: {pgs_file}")

    variants = []
    opener = gzip.open if pgs_file.endswith('.gz') else open

    with opener(pgs_file, 'rt') as f:
        header = None
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            if header is None:
                header = fields
                # Find column indices - prefer harmonized columns
                chr_idx = None
                pos_idx = None

                if 'hm_chr' in header:
                    chr_idx = header.index('hm_chr')
                elif 'chr_name' in header:
                    chr_idx = header.index('chr_name')

                if 'hm_pos' in header:
                    pos_idx = header.index('hm_pos')
                elif 'chr_position' in header:
                    pos_idx = header.index('chr_position')

                try:
                    effect_idx = header.index('effect_allele')
                    weight_idx = header.index('effect_weight')
                    other_idx = header.index('other_allele') if 'other_allele' in header else None
                except ValueError as e:
                    print(f"Error finding required columns: {e}")
                    sys.exit(1)
                continue

            # Skip if no position
            if chr_idx is None or pos_idx is None:
                continue
            if len(fields) <= pos_idx or not fields[pos_idx]:
                continue

            chrom = fields[chr_idx].replace('chr', '')
            pos = int(fields[pos_idx])
            effect = fields[effect_idx]
            weight = float(fields[weight_idx])
            other = fields[other_idx] if other_idx and len(fields) > other_idx and fields[other_idx] else None

            variants.append({
                'chr': chrom,
                'pos': pos,
                'effect_allele': effect,
                'other_allele': other,
                'weight': weight
            })

    print(f"Total variants in PGS file: {len(variants):,}")
    return variants


def load_reference_genome(ref_file, chromosomes_needed):
    """Load reference genome for lookup"""
    print(f"Loading reference genome: {ref_file}")
    print(f"Chromosomes needed: {sorted(set(chromosomes_needed))}")

    ref_dict = {}
    opener = gzip.open if ref_file.endswith('.gz') else open

    current_chr = None
    current_seq = []

    with opener(ref_file, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                # Save previous chromosome
                if current_chr and current_chr in chromosomes_needed:
                    ref_dict[current_chr] = ''.join(current_seq)
                    print(f"  Loaded chr{current_chr}: {len(ref_dict[current_chr]):,} bp")

                # Parse chromosome name
                chr_line = line.strip().split()[0][1:]
                if chr_line.startswith('chr'):
                    current_chr = chr_line[3:]
                else:
                    current_chr = chr_line

                current_seq = []
                if current_chr not in chromosomes_needed:
                    current_chr = None

            elif current_chr:
                current_seq.append(line.strip().upper())

        # Save last chromosome
        if current_chr and current_chr in chromosomes_needed:
            ref_dict[current_chr] = ''.join(current_seq)
            print(f"  Loaded chr{current_chr}: {len(ref_dict[current_chr]):,} bp")

    print(f"Loaded {len(ref_dict)} chromosome sequences")
    return ref_dict


def get_reference_allele(ref_dict, chrom, pos):
    """Get reference allele at position (1-based)"""
    if chrom not in ref_dict:
        return None

    idx = pos - 1
    if idx < 0 or idx >= len(ref_dict[chrom]):
        return None

    return ref_dict[chrom][idx]


def parse_genotype(gt_string):
    """Parse genotype string to get allele indices"""
    if '|' in gt_string:
        alleles = gt_string.split('|')
    elif '/' in gt_string:
        alleles = gt_string.split('/')
    else:
        return None

    try:
        return [int(a) for a in alleles if a != '.']
    except:
        return None


def calculate_pgs_both_methods(vcf_file, pgs_variants, ref_dict, output_prefix, sample_name):
    """
    Calculate PGS using two methods:
    1. ADJUSTED: Includes reference genome lookups (correct for WGS)
    2. UNADJUSTED: VCF data only (for comparison with reference panel)

    Also saves detailed variant-level results
    """

    print(f"\nProcessing VCF: {vcf_file}")

    # Statistics tracking
    stats = {
        'vcf_found': 0,
        'vcf_missing': 0,
        'ref_lookup_success': 0,
        'ref_lookup_failed': 0,
        'allele_mismatch': 0,
        'effect_alt_vcf': 0,
        'effect_ref_vcf': 0,
        'effect_ref_lookup': 0,
        'effect_alt_lookup': 0
    }

    adjusted_score = 0.0
    unadjusted_score = 0.0

    # Store detailed results for output
    detailed_results = []

    # Build lookup dictionary from VCF
    print("Loading VCF variants...")
    vcf_dict = {}
    opener = gzip.open if vcf_file.endswith('.gz') else open

    with opener(vcf_file, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            fields = line.strip().split('\t')

            # Skip malformed lines
            if len(fields) < 10:
                continue

            chrom = fields[0].replace('chr', '')
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            gt_string = fields[9].split(':')[0]

            gt = parse_genotype(gt_string)
            if gt is None:
                continue

            key = (chrom, pos)
            vcf_dict[key] = {
                'ref': ref,
                'alt': alt,
                'gt': gt
            }

    print(f"Loaded {len(vcf_dict):,} variants from VCF")

    # Process each PGS variant
    print("Calculating scores...")
    for variant in pgs_variants:
        chrom = variant['chr']
        pos = variant['pos']
        effect_allele = variant['effect_allele']
        weight = variant['weight']

        key = (chrom, pos)

        if key in vcf_dict:
            # Variant found in VCF
            stats['vcf_found'] += 1

            vcf_info = vcf_dict[key]
            ref = vcf_info['ref']
            alt = vcf_info['alt']
            gt = vcf_info['gt']

            # Check if effect allele matches
            if effect_allele not in [ref, alt]:
                stats['allele_mismatch'] += 1
                continue

            # Calculate dosage
            if effect_allele == alt:
                dosage = sum(gt)
                effect_is_ref = False
                stats['effect_alt_vcf'] += 1
            else:  # effect_allele == ref
                dosage = 2 - sum(gt)
                effect_is_ref = True
                stats['effect_ref_vcf'] += 1

            contribution = dosage * weight

            # Both methods include VCF variants
            adjusted_score += contribution
            unadjusted_score += contribution

            # Save to detailed results
            detailed_results.append({
                'chr': chrom,
                'pos': pos,
                'ref': ref,
                'alt': alt,
                'genotype': f"{gt[0]}/{gt[1]}" if len(gt) >= 2 else str(gt),
                'effect_allele': effect_allele,
                'source': 'vcf',
                'effect_is_ref': effect_is_ref,
                'dosage': dosage,
                'weight': weight,
                'contribution': contribution
            })

        else:
            # Variant NOT in VCF - look up in reference
            stats['vcf_missing'] += 1

            if ref_dict is None:
                stats['ref_lookup_failed'] += 1
                continue

            ref_allele = get_reference_allele(ref_dict, chrom, pos)

            if ref_allele is None:
                stats['ref_lookup_failed'] += 1
                continue

            stats['ref_lookup_success'] += 1

            # Determine dosage based on whether effect allele is reference
            if effect_allele == ref_allele:
                # Homozygous reference, effect allele present at dosage 2
                dosage = 2
                effect_is_ref = True
                stats['effect_ref_lookup'] += 1
            else:
                # Effect allele is NOT reference, dosage = 0
                dosage = 0
                effect_is_ref = False
                stats['effect_alt_lookup'] += 1

            contribution = dosage * weight

            # Only adjusted score includes reference lookups
            adjusted_score += contribution

            # Save to detailed results
            detailed_results.append({
                'chr': chrom,
                'pos': pos,
                'ref': ref_allele,
                'alt': 'N/A',
                'genotype': f"{ref_allele}/{ref_allele}",
                'effect_allele': effect_allele,
                'source': 'ref_lookup',
                'effect_is_ref': effect_is_ref,
                'dosage': dosage,
                'weight': weight,
                'contribution': contribution
            })

    # Save detailed results to file
    detail_file = f"{output_prefix}_detailed_results_v33.txt"
    print(f"\nSaving detailed variant results to {detail_file}...")

    with open(detail_file, 'w') as f:
        f.write("chr\tpos\tref\talt\tgenotype\teffect_allele\tsource\teffect_is_ref\tdosage\tweight\tcontribution\n")
        for result in detailed_results:
            f.write(f"{result['chr']}\t{result['pos']}\t{result['ref']}\t{result['alt']}\t"
                   f"{result['genotype']}\t{result['effect_allele']}\t{result['source']}\t"
                   f"{result['effect_is_ref']}\t{result['dosage']}\t{result['weight']:.6f}\t"
                   f"{result['contribution']:.6f}\n")

    print(f"Saved {len(detailed_results):,} variant results")

    return adjusted_score, unadjusted_score, stats


def extract_pgs_id(pgs_file):
    """Extract PGS ID from filename (e.g., PGS002308 from PGS002308_hmPOS_GRCh38.txt.gz)"""
    basename = os.path.basename(pgs_file)
    match = re.match(r'(PGS\d+)', basename)
    if match:
        return match.group(1)
    # Fallback: use first part before underscore
    return basename.split('_')[0]


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="Calculate Polygenic Risk Scores with both adjusted and unadjusted methods",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full mode with reference genome lookup (recommended for WGS)
  python pgs_calculator.py \\
      --vcf sample.vcf.gz \\
      --pgs-file scores/PGS002308_hmPOS_GRCh38.txt.gz \\
      --reference reference/Homo_sapiens_assembly38.fasta \\
      --build GRCh38 \\
      --output results/PGS002308

  # VCF-only mode (matches pgsc_calc behavior)
  python pgs_calculator.py \\
      --vcf sample.vcf.gz \\
      --pgs-file scores/PGS002308_hmPOS_GRCh38.txt.gz \\
      --no-reference \\
      --output results/PGS002308_vcf_only

Output files:
  {output}_summary_v33.txt       - Tab-separated metrics for notebook
  {output}_detailed_results_v33.txt - Per-variant scoring details
        """,
    )

    parser.add_argument(
        '-v', '--vcf',
        required=True,
        help='Path to VCF file containing your variants'
    )
    parser.add_argument(
        '-p', '--pgs-file',
        required=True,
        help='Path to PGS scoring file (harmonized format)'
    )
    parser.add_argument(
        '-r', '--reference',
        help='Path to reference FASTA file (optional if --no-reference)'
    )
    parser.add_argument(
        '--no-reference',
        action='store_true',
        help='VCF-only mode: skip reference genome lookup (matches pgsc_calc behavior)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output prefix for results (default: derived from PGS ID and sample)'
    )
    parser.add_argument(
        '--build',
        choices=['GRCh37', 'GRCh38'],
        help='Genome build (for documentation purposes)'
    )

    args = parser.parse_args()

    # Validate arguments
    if not args.no_reference and not args.reference:
        parser.error("Either --reference or --no-reference must be specified")

    if args.reference and not os.path.exists(args.reference):
        parser.error(f"Reference file not found: {args.reference}")

    if not os.path.exists(args.vcf):
        parser.error(f"VCF file not found: {args.vcf}")

    if not os.path.exists(args.pgs_file):
        parser.error(f"PGS file not found: {args.pgs_file}")

    return args


def main():
    args = parse_args()

    # Extract PGS ID from filename
    pgs_id = extract_pgs_id(args.pgs_file)

    # Get sample name from VCF
    sample_name = parse_vcf_header(args.vcf)

    # Determine output prefix
    if args.output:
        output_prefix = args.output
    else:
        output_prefix = f"{pgs_id}_{sample_name}"

    # Ensure output directory exists
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    print("="*70)
    print(f"PROPER PGS CALCULATOR v33 - {pgs_id}")
    print("Calculates BOTH adjusted and unadjusted scores")
    print("="*70)

    print(f"\nSample: {sample_name}")
    print(f"VCF: {args.vcf}")
    print(f"PGS File: {args.pgs_file}")
    if args.build:
        print(f"Build: {args.build}")

    if args.no_reference:
        print("\nMode: VCF-only (no reference lookup)")
        print("This matches pgsc_calc behavior")
    else:
        print(f"\nReference: {args.reference}")
        print("Mode: Full (with reference lookup)")

    # Load PGS catalog
    pgs_variants = load_pgs_catalog(args.pgs_file)

    # Get chromosomes needed
    chromosomes = set(v['chr'] for v in pgs_variants)

    # Load reference genome if not in VCF-only mode
    if args.no_reference:
        ref_dict = None
    else:
        ref_dict = load_reference_genome(args.reference, chromosomes)

    # Calculate scores using both methods
    adjusted_score, unadjusted_score, stats = calculate_pgs_both_methods(
        args.vcf, pgs_variants, ref_dict, output_prefix, sample_name
    )

    # Calculate coverage
    total_variants = len(pgs_variants)
    adjusted_coverage = ((stats['vcf_found'] + stats['ref_lookup_success']) / total_variants) * 100
    unadjusted_coverage = (stats['vcf_found'] / total_variants) * 100

    # Print statistics
    print("\n" + "="*70)
    print("SCORING STATISTICS")
    print("="*70)
    print(f"Total variants in PGS file: {total_variants:,}")
    print(f"\nVariant sources:")
    print(f"  Found in VCF: {stats['vcf_found']:,}")
    print(f"  Missing from VCF: {stats['vcf_missing']:,}")
    print(f"    - Looked up in reference: {stats['ref_lookup_success']:,}")
    print(f"    - Reference lookup failed: {stats['ref_lookup_failed']:,}")
    print(f"  Allele mismatches (excluded): {stats['allele_mismatch']:,}")

    print(f"\nEffect allele breakdown:")
    print(f"  From VCF:")
    print(f"    Effect = ALT: {stats['effect_alt_vcf']:,}")
    print(f"    Effect = REF: {stats['effect_ref_vcf']:,}")
    print(f"  From reference lookup:")
    print(f"    Effect = REF (dosage=2): {stats['effect_ref_lookup']:,}")
    print(f"    Effect = ALT (dosage=0): {stats['effect_alt_lookup']:,}")

    # Print scores
    print("\n" + "="*70)
    print("METHOD 1: ADJUSTED SCORE (with reference lookup)")
    print("="*70)
    print(f"Score: {adjusted_score:.4f}")
    print(f"Variants scored: {stats['vcf_found'] + stats['ref_lookup_success']:,}")
    print(f"Coverage: {adjusted_coverage:.1f}%")
    if not args.no_reference:
        print("\nThis is your CORRECT score for WGS data")
        print(f"Includes {stats['ref_lookup_success']:,} homozygous reference sites")

    print("\n" + "="*70)
    print("METHOD 2: UNADJUSTED SCORE (VCF only)")
    print("="*70)
    print(f"Score: {unadjusted_score:.4f}")
    print(f"Variants scored: {stats['vcf_found']:,}")
    print(f"Coverage: {unadjusted_coverage:.1f}%")
    print("\nUse this for comparison with reference panel")
    if not args.no_reference:
        print(f"Excludes {stats['ref_lookup_success']:,} homozygous reference sites")

    difference = adjusted_score - unadjusted_score
    rel_diff = (difference / abs(unadjusted_score) * 100) if unadjusted_score != 0 else 0

    print("\n" + "="*70)
    print("COMPARISON")
    print("="*70)
    print(f"Difference (Adjusted - Unadjusted): {difference:+.4f}")
    print(f"Relative difference: {rel_diff:+.2f}%")

    if not args.no_reference:
        homref_contribution = difference
        if homref_contribution > 0:
            print(f"\nHomozygous reference sites contribute: +{homref_contribution:.4f}")
            print("Impact: Standard methods UNDERESTIMATE your score")
        elif homref_contribution < 0:
            print(f"\nHomozygous reference sites contribute: {homref_contribution:.4f}")
            print("Impact: Standard methods OVERESTIMATE your score")
        else:
            print(f"\nHomozygous reference sites contribute: 0")
            print("Impact: No net effect (contributions cancel)")

    # Save summary results in notebook-compatible format
    summary_file = f"{output_prefix}_summary_v33.txt"
    with open(summary_file, 'w') as f:
        f.write(f"# PGS Calculator v33 Results - {pgs_id}\n")
        f.write(f"# Sample: {sample_name}\n")
        f.write(f"# PGS File: {args.pgs_file}\n")
        f.write(f"# VCF: {args.vcf}\n")
        if args.build:
            f.write(f"# Build: {args.build}\n")
        f.write(f"# Mode: {'VCF-only' if args.no_reference else 'With reference lookup'}\n")
        f.write("\nMetric\tValue\n")
        f.write(f"adjusted_score\t{adjusted_score:.6f}\n")
        f.write(f"unadjusted_score\t{unadjusted_score:.6f}\n")
        f.write(f"difference\t{difference:.6f}\n")
        f.write(f"total_variants\t{total_variants}\n")
        f.write(f"vcf_found\t{stats['vcf_found']}\n")
        f.write(f"ref_lookup_success\t{stats['ref_lookup_success']}\n")
        f.write(f"adjusted_coverage\t{adjusted_coverage:.2f}\n")
        f.write(f"unadjusted_coverage\t{unadjusted_coverage:.2f}\n")
        f.write(f"effect_alt_vcf\t{stats['effect_alt_vcf']}\n")
        f.write(f"effect_ref_vcf\t{stats['effect_ref_vcf']}\n")
        f.write(f"effect_ref_lookup\t{stats['effect_ref_lookup']}\n")
        f.write(f"effect_alt_lookup\t{stats['effect_alt_lookup']}\n")
        f.write(f"allele_mismatch\t{stats['allele_mismatch']}\n")

    print(f"\nResults saved to: {summary_file}")
    print("="*70)


if __name__ == "__main__":
    main()
