#!/usr/bin/env python3
"""
Score mutation size perturbations with Evo2 model.

Tests Evo2's sensitivity to mutation size by scoring mutations with 
different stop codon pattern lengths (1-5 stop codons per insertion).

Questions: 
- Can Evo2 detect minimal perturbations (single TAA) as deleterious?
- Does prediction accuracy improve or degrade with mutation size?

Outputs delta_ll scores for all five size levels.
"""

import pandas as pd
import numpy as np
import torch
from evo2 import Evo2
from tqdm import tqdm
import sys
import gc
import warnings
warnings.filterwarnings('ignore')


def score_sequence(sequence: str, model, device) -> float:
    """Score a single sequence using Evo2 model."""
    try:
        if len(sequence) == 0 or pd.isna(sequence):
            return np.nan
            
        tokens = model.tokenizer.tokenize(sequence)
        
        if not tokens or len(tokens) == 0:
            return np.nan
        
        input_ids = torch.tensor(tokens, dtype=torch.long).unsqueeze(0).to(device)
        
        with torch.no_grad():
            output = model(input_ids)
            if isinstance(output, tuple):
                logits = output[0]
            else:
                logits = output
        
        if isinstance(logits, tuple):
            logits = logits[0]
        
        logprobs = torch.log_softmax(logits, dim=-1)
        logprobs = logprobs[:, :-1]
        input_ids_shifted = input_ids[:, 1:]
        
        seq_logprobs = torch.gather(
            logprobs, 2, input_ids_shifted.unsqueeze(-1)
        ).squeeze(-1)
        
        mean_logprob = seq_logprobs.mean().item()
        return mean_logprob
    
    except Exception as e:
        print(f"Warning - Error scoring sequence: {e}")
        return np.nan


def score_mutation_size_variants(perturbed_csv: str, output_file: str, device: str = 'cuda:0'):
    """
    Score mutations with varying stop codon pattern sizes (1-5 stops per insertion).
    
    Args:
        perturbed_csv: Path to perturbed sequences CSV
        output_file: Output CSV file name
        device: Device to use (cuda:0, cpu, etc.)
    """
    
    print("=" * 80)
    print("EVO2 MUTATION SIZE PERTURBATION SCORING")
    print("=" * 80)
    
    # Load perturbed sequences
    print(f"\n[1/4] Loading mutation size perturbations...")
    print(f"      Input: {perturbed_csv}")
    genes_df = pd.read_csv(perturbed_csv)
    print(f"      Loaded {len(genes_df)} genes")
    
    # Load Evo2 model
    print(f"\n[2/4] Loading Evo2 model...")
    device = torch.device(device if torch.cuda.is_available() else 'cpu')
    print(f"      Device: {device}")
    
    model = Evo2('evo2_7b')
    print(f"      Model loaded successfully")
    
    # Process each gene
    print(f"\n[3/4] Scoring sequences...")
    results = []
    
    for idx, row in tqdm(genes_df.iterrows(), total=len(genes_df), desc="Progress"):
        gene_id = row['gene_id']
        wt_sequence = row['wt_sequence']
        
        try:
            # Score wild-type
            wt_score = score_sequence(wt_sequence, model, device)
            
            # Score mutations of different sizes
            size1_score = score_sequence(row['mut_size1_sequence'], model, device) if pd.notna(row.get('mut_size1_sequence')) else np.nan
            size2_score = score_sequence(row['mut_size2_sequence'], model, device) if pd.notna(row.get('mut_size2_sequence')) else np.nan
            size3_score = score_sequence(row['mut_size3_sequence'], model, device) if pd.notna(row.get('mut_size3_sequence')) else np.nan
            size4_score = score_sequence(row['mut_size4_sequence'], model, device) if pd.notna(row.get('mut_size4_sequence')) else np.nan
            size5_score = score_sequence(row['mut_size5_sequence'], model, device) if pd.notna(row.get('mut_size5_sequence')) else np.nan
            
            # Calculate delta log-likelihood
            delta_ll_size1 = size1_score - wt_score if not np.isnan(size1_score) and not np.isnan(wt_score) else np.nan
            delta_ll_size2 = size2_score - wt_score if not np.isnan(size2_score) and not np.isnan(wt_score) else np.nan
            delta_ll_size3 = size3_score - wt_score if not np.isnan(size3_score) and not np.isnan(wt_score) else np.nan
            delta_ll_size4 = size4_score - wt_score if not np.isnan(size4_score) and not np.isnan(wt_score) else np.nan
            delta_ll_size5 = size5_score - wt_score if not np.isnan(size5_score) and not np.isnan(wt_score) else np.nan
            
            results.append({
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
                'wt_logprob': wt_score,
                'mut_size1_logprob': size1_score,
                'mut_size2_logprob': size2_score,
                'mut_size3_logprob': size3_score,
                'mut_size4_logprob': size4_score,
                'mut_size5_logprob': size5_score,
                'delta_ll_size1': delta_ll_size1,
                'delta_ll_size2': delta_ll_size2,
                'delta_ll_size3': delta_ll_size3,
                'delta_ll_size4': delta_ll_size4,
                'delta_ll_size5': delta_ll_size5,
                'stop_pattern_size1': row.get('stop_pattern_size1', 'TAA'),
                'stop_pattern_size2': row.get('stop_pattern_size2', 'TAATAA'),
                'stop_pattern_size3': row.get('stop_pattern_size3', 'TAATAATAA'),
                'stop_pattern_size4': row.get('stop_pattern_size4', 'TAATAATAATAG'),
                'stop_pattern_size5': row.get('stop_pattern_size5', 'TAATAATAATAGTGA'),
                'inserted_percent_cds_size1': row.get('inserted_percent_cds_size1', np.nan),
                'inserted_percent_cds_size2': row.get('inserted_percent_cds_size2', np.nan),
                'inserted_percent_cds_size3': row.get('inserted_percent_cds_size3', np.nan),
                'inserted_percent_cds_size4': row.get('inserted_percent_cds_size4', np.nan),
                'inserted_percent_cds_size5': row.get('inserted_percent_cds_size5', np.nan),
            })
        
        except Exception as e:
            print(f"      Warning - Error processing {gene_id}: {e}")
            results.append({
                'gene_id': gene_id,
                'locus_tag': row.get('locus_tag', ''),
                'gene_name': row.get('gene_name', ''),
                'cds_length': row.get('cds_length', np.nan),
                'wt_logprob': np.nan,
                'mut_size1_logprob': np.nan,
                'mut_size2_logprob': np.nan,
                'mut_size3_logprob': np.nan,
                'mut_size4_logprob': np.nan,
                'mut_size5_logprob': np.nan,
                'delta_ll_size1': np.nan,
                'delta_ll_size2': np.nan,
                'delta_ll_size3': np.nan,
                'delta_ll_size4': np.nan,
                'delta_ll_size5': np.nan,
                'stop_pattern_size1': 'TAA',
                'stop_pattern_size2': 'TAATAA',
                'stop_pattern_size3': 'TAATAATAA',
                'stop_pattern_size4': 'TAATAATAATAG',
                'stop_pattern_size5': 'TAATAATAATAGTGA',
                'inserted_percent_cds_size1': np.nan,
                'inserted_percent_cds_size2': np.nan,
                'inserted_percent_cds_size3': np.nan,
                'inserted_percent_cds_size4': np.nan,
                'inserted_percent_cds_size5': np.nan,
            })
        
        # Memory cleanup
        if idx % 100 == 0:
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    print(f"\n[4/4] Saving results...")
    output_df.to_csv(output_file, index=False)
    print(f"      Output: {output_file}")
    
    # Print summary statistics
    print(f"\n" + "=" * 80)
    print(f"SUMMARY STATISTICS")
    print(f"=" * 80)
    print(f"Total genes processed: {len(output_df)}")
    
    valid_counts = {
        'size1': output_df['delta_ll_size1'].notna().sum(),
        'size2': output_df['delta_ll_size2'].notna().sum(),
        'size3': output_df['delta_ll_size3'].notna().sum(),
        'size4': output_df['delta_ll_size4'].notna().sum(),
        'size5': output_df['delta_ll_size5'].notna().sum(),
    }
    
    print(f"\nGenes with valid scores:")
    for size, count in valid_counts.items():
        print(f"  {size}: {count}")
    
    if valid_counts['size1'] > 0:
        print(f"\n" + "-" * 80)
        print(f"DELTA LOG-LIKELIHOOD BY MUTATION SIZE")
        print(f"-" * 80)
        
        for i in range(1, 6):
            vals = output_df[f'delta_ll_size{i}'].dropna()
            pattern = output_df[f'stop_pattern_size{i}'].iloc[0] if f'stop_pattern_size{i}' in output_df.columns else f'size{i}'
            print(f"\nSIZE {i} ({pattern}):")
            print(f"  Mean:   {vals.mean():.6f}")
            print(f"  Median: {vals.median():.6f}")
            print(f"  Std:    {vals.std():.6f}")
        
        # Statistical interpretation
        print(f"\n" + "-" * 80)
        print(f"INTERPRETATION")
        print(f"-" * 80)
        
        means = [output_df[f'delta_ll_size{i}'].dropna().mean() for i in range(1, 6)]
        
        print(f"\nMean delta_ll by size:")
        for i, m in enumerate(means, 1):
            print(f"  Size {i}: {m:.6f}")
        
        # Check if smaller mutations still produce detectable signal
        if all(m < 0 for m in means):
            print(f"\n✓ All mutation sizes produce negative delta_ll (detectable as deleterious)")
        else:
            print(f"\n✗ Some mutation sizes do not produce negative delta_ll")
        
        # Correlation with mutation size
        from scipy.stats import spearmanr
        all_delta_ll = pd.concat([output_df[f'delta_ll_size{i}'].dropna() for i in range(1, 6)])
        all_sizes = pd.concat([pd.Series([i] * output_df[f'delta_ll_size{i}'].notna().sum()) for i in range(1, 6)]).reset_index(drop=True)
        
        corr, p_value = spearmanr(all_sizes, all_delta_ll)
        print(f"\nSpearman correlation (pattern size vs delta_ll):")
        print(f"  ρ = {corr:.4f}, p-value = {p_value:.4e}")
    
    print(f"\n" + "=" * 80)
    print(f"SCORING COMPLETE")
    print(f"=" * 80 + "\n")
    
    return output_df


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python score_mutation_size_variants.py <input.csv> [output.csv] [device]")
        print("Example: python score_mutation_size_variants.py genes_mutation_size_perturbed.csv mutation_size_scores.csv cuda:0")
        sys.exit(1)
    
    perturbed_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'essentiality_scores_mutation_size.csv'
    device = sys.argv[3] if len(sys.argv) > 3 else 'cuda:0'
    
    scores_df = score_mutation_size_variants(perturbed_csv, output_file, device=device)