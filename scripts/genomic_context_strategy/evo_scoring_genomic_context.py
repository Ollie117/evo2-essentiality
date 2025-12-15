#!/usr/bin/env python3
"""
Score perturbed and wild-type sequences using Evo2 model - Genomic Context Strategy.

Loads sequences (wild-type and mutant) across 5 context levels and calculates 
log-likelihood scores for each using the Evo2 7B model. Computes delta log-likelihood 
(Δℓ) which serves as a proxy for fitness impact and essentiality prediction.

Context levels:
- Level 1: 512bp each side (1024bp total)
- Level 2: 1024bp each side (2048bp total)
- Level 3: 2048bp each side (4096bp total)
- Level 4: 4096bp each side (8192bp total) - BASELINE
- Level 5: 8192bp each side (16384bp total)
"""

import pandas as pd
import numpy as np
import torch
from evo2 import Evo2
from tqdm import tqdm
import sys

CONTEXT_LEVELS = [1, 2, 3, 4, 5]

def score_sequence(sequence: str, model: Evo2, device: str = 'cuda:0') -> float:
    """
    Score a single sequence using Evo2 model.
    
    Args:
        sequence: DNA sequence to score
        model: Loaded Evo2 model
        device: Device to use (cuda:0, cpu, etc.)
    
    Returns:
        Mean log-likelihood score for the sequence
    """
    try:
        # Tokenize sequence
        tokens = model.tokenizer.tokenize(sequence)
        input_ids = torch.tensor(tokens, dtype=torch.long).unsqueeze(0).to(device)
        
        # Get logits
        with torch.no_grad():
            output = model(input_ids)
            # Handle tuple output - logits is typically first element
            if isinstance(output, tuple):
                logits = output[0]
            else:
                logits = output
        
        # Convert to log probabilities
        logprobs = torch.log_softmax(logits, dim=-1)
        
        # Get log probability of actual tokens
        # Remove last position (lookahead) and first token (BOS)
        logprobs = logprobs[:, :-1]
        input_ids = input_ids[:, 1:]
        
        # Gather log probabilities for actual tokens
        seq_logprobs = torch.gather(
            logprobs,
            2,
            input_ids.unsqueeze(-1)
        ).squeeze(-1)
        
        # Calculate mean log-likelihood
        mean_logprob = seq_logprobs.mean().item()
        return mean_logprob
    
    except Exception as e:
        print(f"Error scoring sequence: {e}")
        return np.nan

def score_perturbed_sequences(perturbed_csv: str, output_file: str = 'essentiality_scores_genomic_context.csv',
                              batch_size: int = 1, device: str = 'cuda:0'):
    """
    Score wild-type and mutant sequences across all context levels, calculate delta log-likelihood.
    
    Args:
        perturbed_csv: Path to perturbed sequences CSV (from generate_genomic_context_sequences.py)
        output_file: Output CSV file name
        batch_size: Batch size for processing (currently processes one at a time)
        device: Device to use (cuda:0, cpu, etc.)
    """
    
    # Load perturbed sequences
    print(f"Loading perturbed sequences from: {perturbed_csv}")
    genes_df = pd.read_csv(perturbed_csv)
    print(f"Loaded {len(genes_df)} genes")
    
    # Load Evo2 model
    print("\nLoading Evo2 model (7B variant)...")
    device = torch.device(device if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")
    
    model = Evo2('evo2_7b')
    print("Model loaded successfully\n")
    
    # Process each gene
    results = []
    
    for idx, row in tqdm(genes_df.iterrows(), total=len(genes_df), desc="Scoring sequences"):
        gene_id = row['gene_id']
        
        try:
            result = {
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
            }
            
            # Score each context level
            for level in CONTEXT_LEVELS:
                wt_sequence = row[f'wt_context{level}_sequence']
                mut_sequence = row[f'mut_context{level}_sequence']
                
                # Score wild-type
                wt_score = score_sequence(wt_sequence, model, device=str(device))
                
                # Score mutant
                mut_score = score_sequence(mut_sequence, model, device=str(device))
                
                # Calculate delta log-likelihood
                delta_ll = mut_score - wt_score
                
                result[f'wt_context{level}_logprob'] = wt_score
                result[f'mut_context{level}_logprob'] = mut_score
                result[f'delta_ll_context{level}'] = delta_ll
                result[f'context{level}_bp'] = row[f'context{level}_bp']
            
            results.append(result)
        
        except Exception as e:
            print(f"Error processing {gene_id}: {e}")
            result = {
                'gene_id': gene_id,
                'locus_tag': row['locus_tag'],
                'gene_name': row['gene_name'],
                'cds_length': row['cds_length'],
            }
            for level in CONTEXT_LEVELS:
                result[f'wt_context{level}_logprob'] = np.nan
                result[f'mut_context{level}_logprob'] = np.nan
                result[f'delta_ll_context{level}'] = np.nan
                result[f'context{level}_bp'] = np.nan
            results.append(result)
    
    # Create output dataframe
    output_df = pd.DataFrame(results)
    
    # Save results
    output_df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")
    
    # Print summary statistics
    print(f"\nSummary Statistics:")
    print(f"  Total genes processed: {len(output_df)}")
    
    for level in CONTEXT_LEVELS:
        col = f'delta_ll_context{level}'
        valid_scores = output_df[col].notna().sum()
        if valid_scores > 0:
            print(f"\n  Context Level {level} ({512 * (2 ** (level-1))}bp each side):")
            print(f"    Genes with valid scores: {valid_scores}")
            print(f"    Mean Δℓ: {output_df[col].mean():.6f}")
            print(f"    Median Δℓ: {output_df[col].median():.6f}")
            print(f"    Negative Δℓ (deleterious): {(output_df[col] < 0).sum()} genes")
            print(f"    Positive Δℓ (tolerated): {(output_df[col] > 0).sum()} genes")
    
    # Show example predictions
    print(f"\nExample predictions (first 5 genes, context level 4 - baseline):")
    print(output_df.head(5)[['gene_id', 'gene_name', 'wt_context4_logprob', 'mut_context4_logprob', 'delta_ll_context4']].to_string())
    
    return output_df

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python evo_scoring_genomic_context.py <perturbed_sequences.csv> [output_file] [device]")
        print("Example: python evo_scoring_genomic_context.py genes_genomic_context.csv essentiality_scores_genomic_context.csv cuda:0")
        sys.exit(1)
    
    perturbed_csv = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'essentiality_scores_genomic_context.csv'
    device = sys.argv[3] if len(sys.argv) > 3 else 'cuda:0'
    
    scores_df = score_perturbed_sequences(perturbed_csv, output_file, device=device)