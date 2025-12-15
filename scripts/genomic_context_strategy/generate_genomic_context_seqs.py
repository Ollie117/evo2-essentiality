#!/usr/bin/env python3
"""
Generate perturbed sequences with varying genomic context windows.

Tests whether Evo2 relies on distant genomic elements by varying the amount
of flanking sequence provided to the model.

Context levels:
- Level 1: 512bp each side (1024bp total)
- Level 2: 1024bp each side (2048bp total)
- Level 3: 2048bp each side (4096bp total)
- Level 4: 4096bp each side (8192bp total) - BASELINE
- Level 5: 8192bp each side (16384bp total)
"""

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import re
import os

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
GENOME_FILE = "data/raw/GCF_000195955.2_ASM19595v2_genomic.fna"
ANNOTATION_FILE = "data/raw/genomic.gtf"
TNSEQ_FILE = "data/processed/tnseq_clean.csv"

# Output file
OUTPUT_FILE = "data/processed/genes_genomic_context.csv"

# Context levels (bp on each side)
CONTEXT_LEVELS = {
    1: 512,    # 1024bp total
    2: 1024,   # 2048bp total
    3: 2048,   # 4096bp total
    4: 4096,   # 8192bp total - BASELINE
    5: 8192,   # 16384bp total
}

# Perturbation parameters (baseline strategy)
STOP_PATTERN = "TAATAATAATAGTGA"
INSERTION_INTERVAL = 12

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading genome sequence...")
genome_record = next(SeqIO.parse(GENOME_FILE, "fasta"))
genome_seq = str(genome_record.seq)
genome_length = len(genome_seq)
print(f"  Genome length: {genome_length:,} bp")

print("\nLoading gene annotations from GTF...")
genes = []
with open(ANNOTATION_FILE, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if len(fields) >= 9 and fields[2] == 'gene':
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            
            # Parse GTF attributes (format: key "value"; key "value";)
            attrs_str = fields[8]
            
            # Extract gene_id
            gene_id_match = re.search(r'gene_id "([^"]+)"', attrs_str)
            gene_id = gene_id_match.group(1) if gene_id_match else ''
            
            # Extract locus_tag
            locus_tag_match = re.search(r'locus_tag "([^"]+)"', attrs_str)
            locus_tag = locus_tag_match.group(1) if locus_tag_match else gene_id
            
            # Extract gene name
            gene_name_match = re.search(r'gene "([^"]+)"', attrs_str)
            gene_name = gene_name_match.group(1) if gene_name_match else ''
            
            if gene_id:
                genes.append({
                    'gene_id': gene_id,
                    'locus_tag': locus_tag,
                    'gene_name': gene_name,
                    'start': start,
                    'end': end,
                    'strand': strand
                })

genes_df = pd.DataFrame(genes)
print(f"  Loaded {len(genes_df)} genes")

print("\nLoading TnSeq classifications...")
tnseq_df = pd.read_csv(TNSEQ_FILE)
print(f"  Loaded {len(tnseq_df)} classified genes")

# Filter to genes with essentiality classifications
genes_df = genes_df[genes_df['gene_id'].isin(tnseq_df['gene_id'])]
print(f"  Genes with essentiality data: {len(genes_df)}")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_sequence_with_context(start, end, strand, context_bp):
    """Extract gene sequence with specified flanking context."""
    # Calculate context boundaries
    context_start = max(0, start - 1 - context_bp)
    context_end = min(genome_length, end + context_bp)
    
    # Extract sequence
    seq = genome_seq[context_start:context_end]
    
    # Reverse complement if on negative strand
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    
    # Calculate CDS position within the extracted sequence
    if strand == '+':
        cds_start_in_seq = (start - 1) - context_start
        cds_end_in_seq = end - context_start
    else:
        cds_start_in_seq = context_end - end
        cds_end_in_seq = context_end - (start - 1)
    
    return seq, cds_start_in_seq, cds_end_in_seq, context_end - context_start

def insert_stop_codons(sequence, cds_start, cds_end, pattern, interval):
    """Insert stop codon pattern throughout CDS at specified interval."""
    cds_seq = sequence[cds_start:cds_end]
    cds_length = len(cds_seq)
    
    # Calculate insertion positions
    insertions = []
    pos = interval
    while pos < cds_length:
        insertions.append(pos)
        pos += interval
    
    # Insert patterns (from end to start to preserve positions)
    mutant_cds = cds_seq
    for pos in reversed(insertions):
        mutant_cds = mutant_cds[:pos] + pattern + mutant_cds[pos:]
    
    # Reconstruct full sequence
    mutant_seq = sequence[:cds_start] + mutant_cds + sequence[cds_end:]
    
    return mutant_seq, len(insertions)

# ============================================================================
# GENERATE SEQUENCES
# ============================================================================

print("\nGenerating sequences with varying context...")
results = []

for idx, gene in genes_df.iterrows():
    gene_id = gene['gene_id']
    start = gene['start']
    end = gene['end']
    strand = gene['strand']
    cds_length = end - start + 1
    
    gene_result = {
        'gene_id': gene_id,
        'locus_tag': gene['locus_tag'],
        'gene_name': gene['gene_name'],
        'cds_length': cds_length,
    }
    
    # Generate sequences for each context level
    for level, context_bp in CONTEXT_LEVELS.items():
        # Get wild-type sequence with context
        wt_seq, cds_start, cds_end, total_length = get_sequence_with_context(
            start, end, strand, context_bp
        )
        
        # Generate mutant sequence
        mut_seq, num_insertions = insert_stop_codons(
            wt_seq, cds_start, cds_end, STOP_PATTERN, INSERTION_INTERVAL
        )
        
        # Calculate metrics
        inserted_bp = num_insertions * len(STOP_PATTERN)
        inserted_percent = (inserted_bp / cds_length) * 100
        
        # Store results
        gene_result[f'wt_context{level}_sequence'] = wt_seq
        gene_result[f'mut_context{level}_sequence'] = mut_seq
        gene_result[f'context{level}_bp'] = context_bp
        gene_result[f'context{level}_total_bp'] = total_length
        gene_result[f'num_insertions_context{level}'] = num_insertions
    
    # Common metadata
    gene_result['stop_pattern'] = STOP_PATTERN
    gene_result['insertion_interval'] = INSERTION_INTERVAL
    gene_result['perturbation_method'] = 'genomic_context'
    
    results.append(gene_result)
    
    if (idx + 1) % 500 == 0:
        print(f"  Processed {idx + 1}/{len(genes_df)} genes...")

# ============================================================================
# SAVE RESULTS
# ============================================================================

print(f"\nSaving results to {OUTPUT_FILE}...")
results_df = pd.DataFrame(results)
results_df.to_csv(OUTPUT_FILE, index=False)

print(f"\nGenerated sequences for {len(results_df)} genes")
print(f"Context levels tested:")
for level, context_bp in CONTEXT_LEVELS.items():
    print(f"  Level {level}: {context_bp}bp each side ({context_bp * 2}bp total)")

print("\n=== Sequence Generation Complete ===")