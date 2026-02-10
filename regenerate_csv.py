#!/usr/bin/env python3
"""Regenerate comprehensive_bias_analysis.csv with both Trent and Rowen data."""

import pandas as pd
from pathlib import Path
import re

# PGS ID to trait mapping
PGS_TRAITS = {
    'PGS002308': 'Type 2 Diabetes',
    'PGS004034': "Alzheimer's Disease",
    'PGS000027': 'BMI',
    'PGS004237': 'CAD',
    'PGS004696': 'CHD'
}

def parse_summary_file(filepath: Path) -> dict:
    """Parse a v33 summary file and extract metrics."""
    metrics = {}
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            if '\t' in line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    key, value = parts[0], parts[1]
                    try:
                        metrics[key] = float(value)
                    except ValueError:
                        metrics[key] = value
    return metrics

def calculate_bias(adjusted: float, unadjusted: float) -> float:
    """Calculate homozygous reference bias percentage.

    Formula from PRS_Methodology_Experiment.ipynb:
    (adjusted - unadjusted) / abs(unadjusted) * 100

    Using abs(unadjusted) ensures consistent sign interpretation:
    - Positive bias = adjusted > unadjusted (VCF-only underestimates magnitude)
    - Negative bias = adjusted < unadjusted (VCF-only overestimates magnitude)
    """
    if unadjusted == 0:
        return 0.0
    return ((adjusted - unadjusted) / abs(unadjusted)) * 100

def get_interpretation(bias_pct: float) -> str:
    """Get interpretation string based on bias percentage.

    Positive bias = adjusted > unadjusted = VCF-only UNDERESTIMATES
    Negative bias = adjusted < unadjusted = VCF-only OVERESTIMATES
    """
    abs_bias = abs(bias_pct)
    if abs_bias < 5:
        return "✓ Minimal bias"
    elif abs_bias < 50:
        return "⚠️ Underestimate" if bias_pct > 0 else "⚠️ Overestimate"
    else:
        return "⚠️ SEVERE UNDERESTIMATE" if bias_pct > 0 else "⚠️ SEVERE OVERESTIMATE"

def main():
    base_path = Path('results/methodology_comparison/custom')
    samples = {'trent': 'Trent', 'rowan': 'Rowen'}

    results = []

    for sample_dir, sample_name in samples.items():
        sample_path = base_path / sample_dir
        if not sample_path.exists():
            print(f"Warning: {sample_path} does not exist")
            continue

        for pgs_id, trait in PGS_TRAITS.items():
            summary_file = sample_path / f'{pgs_id}_summary_v33.txt'
            if not summary_file.exists():
                print(f"Warning: {summary_file} does not exist")
                continue

            metrics = parse_summary_file(summary_file)

            adjusted = metrics.get('adjusted_score', 0)
            unadjusted = metrics.get('unadjusted_score', 0)
            bias_pct = calculate_bias(adjusted, unadjusted)

            row = {
                'sample': sample_name,
                'pgs_id': pgs_id,
                'trait': trait,
                'adjusted_score': adjusted,
                'unadjusted_score': unadjusted,
                'homref_bias_pct': bias_pct,
                'variants_from_vcf': int(metrics.get('vcf_found', 0)),
                'variants_from_ref': int(metrics.get('ref_lookup_success', 0)),
                'pgsc_raw_score': '',  # pgsc_calc data not available for Rowen
                'pgsc_z_score': '',
                'pgsc_percentile': '',
                'raw_score_diff': '',
                'raw_score_diff_pct': '',
                'interpretation': get_interpretation(bias_pct)
            }
            results.append(row)

    df = pd.DataFrame(results)

    # Sort by sample then trait for consistent ordering
    df = df.sort_values(['sample', 'trait']).reset_index(drop=True)

    output_path = Path('results/methodology_comparison/comprehensive_bias_analysis.csv')
    df.to_csv(output_path, index=False)
    print(f"Generated {output_path} with {len(df)} rows")
    print(f"\nSamples: {df['sample'].unique().tolist()}")
    print(f"Traits: {df['trait'].unique().tolist()}")
    print(f"\nPreview:")
    print(df[['sample', 'trait', 'adjusted_score', 'homref_bias_pct', 'interpretation']].to_string())

if __name__ == '__main__':
    main()
