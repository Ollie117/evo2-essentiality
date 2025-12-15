#!/usr/bin/env python3
"""
Generate mutations with varying stop codon pattern sizes to test model sensitivity.

Creates perturbations where stop codons are inserted at 12bp intervals,
but with different pattern sizes (1, 2, 3, 4, or 5 stop codons per insertion).

Tests: Can Evo2 detect smaller/subtler mutations as deleterious?
Does it require large perturbations or can it recognize minimal disruptions?

Keeps interval constant at 12bp (as per Evo2 paper), varies pattern size.
"""

import pandas as pd
import numpy as np
import sys

# Define stop codon patterns of increasing size
STOP_PATTERNS = {
    'size1': 'TAA',                    # 3bp, 1 stop
    'size2': 'TAATAA',                 # 6bp, 2 stops
    'size3': 'TAATAATAA',              # 9bp, 3 stops
    'size4': 'TAATAATAATAG',           # 12bp, 4 stops
    'size5': 'TAATAATAATAGTGA',        # 15bp, 5 stops (baseline/paper method)
}

def insert_stops_at_intervals(sequence: str, offset: int = 12, stop_pattern: str = 'TAA') -> str:
    """
    Insert stop codon pattern at regular intervals throughout sequence.
    
    Args:
        sequence: DNA sequence (CDS only)
        offset: Interval between insertions (default: 12bp, as per paper)
        stop_pattern: Stop codon pattern to insert
    
    Returns:
        Sequence with stops inserted at intervals
    """
    if len(sequence) < offset:
        return stop_pattern + sequence
    
    mutated = []
    for i in range(0, len(sequence), offset):
        mutated.append(stop_pattern)
        mutated.append(sequence[i:i+offset])
    
    return ''.join(mutated)


def generate_mutation_size_perturbations(genes_csv: str, output_file: str = 'genes_mutation_size_perturbed.csv'):
    """
    Generate mutations with varying stop codon pattern sizes.
    
    For each gene, create five perturbations with stops inserted every 12bp:
    - Size 1: TAA (3bp, 1 stop) - minimal
    - Size 2: TAATAA (6bp, 2 stops)
    - Size 3: TAATAATAA (9bp, 3 stops)
    - Size 4: TAATAATAATAG (12bp, 4 stops)
    - Size 5: TAATAATAATAGTGA (15bp, 5 stops) - baseline/paper method
    
    Tests: How sensitive is Evo2 to mutation size?
    
    Args:
        genes_csv: Path to genes_sequences_8k.csv (with wt_sequence and context info)
        output_file: Output CSV file path
    """
    
    print("=" * 80)
    print("MUTATION SIZE PERTURBATION GENERATION")
    print("=" * 80)
    
    print(f"\nStop codon patterns:")
    for name, pattern in STOP_PATTERNS.items():
        print(f"  {name}: {pattern} ({len(pattern)}bp, {pattern.count('TAA') + pattern.count('TAG') + pattern.count('TGA')} stops)")
    
    # Load gene sequences
    print(f"\n[1/2] Loading gene sequences from: {genes_csv}")
    genes_df = pd.read_csv(genes_csv)
    print(f"      Loaded {len(genes_df)} genes")
    
    # Generate perturbations
    print(f"\n[2/2] Generating mutation size perturbations...")
    results = []
    errors = []
    
    for idx, row in genes_df.iterrows():
        if idx % 500 == 0:
            print(f"      Processing gene {idx}/{len(genes_df)}")
        
        gene_id = row['gene_id']
        wt_sequence = row['sequence']
        
        try:
            cds_start = int(row['cds_start_in_seq'])
            cds_end = int(row['cds_end_in_seq'])
            cds_sequence = wt_sequence[cds_start:cds_end]
            cds_length = len(cds_sequence)
            
            # Extract context
            context_left = wt_sequence[:cds_start]
            context_right = wt_sequence[cds_end:]
            
            # Create mutations with different pattern sizes (all at 12bp intervals)
            mut_size1_cds = insert_stops_at_intervals(cds_sequence, offset=12, stop_pattern=STOP_PATTERNS['size1'])
            mut_size2_cds = insert_stops_at_intervals(cds_sequence, offset=12, stop_pattern=STOP_PATTERNS['size2'])
            mut_size3_cds = insert_stops_at_intervals(cds_sequence, offset=12, stop_pattern=STOP_PATTERNS['size3'])
            mut_size4_cds = insert_stops_at_intervals(cds_sequence, offset=12, stop_pattern=STOP_PATTERNS['size4'])
            mut_size5_cds = insert_stops_at_intervals(cds_sequence, offset=12, stop_pattern=STOP_PATTERNS['size5'])
            
            # Reconstruct full sequences with genomic context
            mut_size1_full = context_left + mut_size1_cds + context_right
            mut_size2_full = context_left + mut_size2_cds + context_right
            mut_size3_full = context_left + mut_size3_cds + context_right
            mut_size4_full = context_left + mut_size4_cds + context_right
            mut_size5_full = context_left + mut_size5_cds + context_right
            
            # Calculate mutation burden (inserted bp as % of original CDS)
            num_intervals = (cds_length // 12) + 1
            
            inserted_bp_1 = num_intervals * len(STOP_PATTERNS['size1'])
            inserted_bp_2 = num_intervals * len(STOP_PATTERNS['size2'])
            inserted_bp_3 = num_intervals * len(STOP_PATTERNS['size3'])
            inserted_bp_4 = num_intervals * len(STOP_PATTERNS['size4'])
            inserted_bp_5 = num_intervals * len(STOP_PATTERNS['size5'])
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'context_bp': row.get('context_bp', 4096),
                'wt_sequence': wt_sequence,
                'mut_size1_sequence': mut_size1_full,
                'mut_size2_sequence': mut_size2_full,
                'mut_size3_sequence': mut_size3_full,
                'mut_size4_sequence': mut_size4_full,
                'mut_size5_sequence': mut_size5_full,
                'stop_pattern_size1': STOP_PATTERNS['size1'],
                'stop_pattern_size2': STOP_PATTERNS['size2'],
                'stop_pattern_size3': STOP_PATTERNS['size3'],
                'stop_pattern_size4': STOP_PATTERNS['size4'],
                'stop_pattern_size5': STOP_PATTERNS['size5'],
                'inserted_bp_size1': inserted_bp_1,
                'inserted_bp_size2': inserted_bp_2,
                'inserted_bp_size3': inserted_bp_3,
                'inserted_bp_size4': inserted_bp_4,
                'inserted_bp_size5': inserted_bp_5,
                'inserted_percent_cds_size1': (inserted_bp_1 / cds_length) * 100,
                'inserted_percent_cds_size2': (inserted_bp_2 / cds_length) * 100,
                'inserted_percent_cds_size3': (inserted_bp_3 / cds_length) * 100,
                'inserted_percent_cds_size4': (inserted_bp_4 / cds_length) * 100,
                'inserted_percent_cds_size5': (inserted_bp_5 / cds_length) * 100,
                'num_stops_size1': 1,
                'num_stops_size2': 2,
                'num_stops_size3': 3,
                'num_stops_size4': 4,
                'num_stops_size5': 5,
                'perturbation_interval': 12,
                'perturbation_method': 'mutation_size'
            })
        
        except Exception as e:
            print(f"      Warning - Error processing {gene_id}: {e}")
            errors.append({'gene_id': gene_id, 'error': str(e)})
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    output_df.to_csv(output_file, index=False)
    print(f"\n      Results saved to: {output_file}")
    
    # Print summary statistics
    print(f"\n" + "=" * 80)
    print(f"SUMMARY STATISTICS")
    print(f"=" * 80)
    print(f"Total genes processed: {len(genes_df)}")
    print(f"Genes successfully perturbed: {len(output_df)}")
    print(f"Errors: {len(errors)}")
    
    print(f"\nMutation burden by size level (average across genes):")
    print(f"  Size 1 (TAA, 1 stop):              {output_df['inserted_percent_cds_size1'].mean():.2f}% of CDS inserted")
    print(f"  Size 2 (TAATAA, 2 stops):          {output_df['inserted_percent_cds_size2'].mean():.2f}% of CDS inserted")
    print(f"  Size 3 (TAATAATAA, 3 stops):       {output_df['inserted_percent_cds_size3'].mean():.2f}% of CDS inserted")
    print(f"  Size 4 (TAATAATAATAG, 4 stops):    {output_df['inserted_percent_cds_size4'].mean():.2f}% of CDS inserted")
    print(f"  Size 5 (TAATAATAATAGTGA, 5 stops): {output_df['inserted_percent_cds_size5'].mean():.2f}% of CDS inserted")
    
    print(f"\nExample genes:")
    for idx, row in output_df.head(2).iterrows():
        print(f"\n  {row['gene_id']} ({row['gene_name']}):")
        print(f"    CDS length: {row['cds_length']} bp")
        print(f"    WT sequence length: {len(row['wt_sequence'])} bp")
        print(f"    Mutant lengths: size1={len(row['mut_size1_sequence'])}, size5={len(row['mut_size5_sequence'])} bp")
    
    print(f"\n" + "=" * 80)
    print(f"GENERATION COMPLETE")
    print(f"=" * 80 + "\n")
    
    return output_df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python generate_mutation_size_perturbations.py <genes_sequences_8k.csv> [output_file]")
        print("Example: python generate_mutation_size_perturbations.py data/processed/genes_sequences_8k.csv genes_mutation_size_perturbed.csv")
        sys.exit(1)
    
    genes_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'genes_mutation_size_perturbed.csv'
    
    perturbed_df = generate_mutation_size_perturbations(genes_csv, output_file)