"""
Extract capsid proteins from Orthoretrovirinae Gag proteins in SwissProt.

Workflow:
1. Download all SwissProt Gag proteins from Orthoretrovirinae
2. Extract capsid features (complete, non-nucleocapsid)
3. Extract sequences and metadata
4. Find unique sequences and aggregate metadata
"""

import os
import json
import sys

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from utils import uniprot_utils
from utils import metadata_utils
from utils import fasta_utils
from utils import mmseqs2_utils


def main():
    """
    Main workflow function to extract and process capsid proteins from Orthoretrovirinae Gag proteins.
    
    This function orchestrates the complete pipeline:
    1. Downloads UniProt SwissProt entries for Orthoretrovirinae Gag proteins
    2. Extracts capsid features (complete, non-nucleocapsid) from each entry
    3. Aggregates metadata by unique sequences
    4. Clusters sequences using MMseqs2
    5. Saves results in multiple formats (JSON, TSV, FASTA)
    
    All output files are saved in the 'outputs' directory.
    """
    # Set up output directories
    output_dir = 'outputs'
    mmseqs2_output_dir = os.path.join(output_dir, 'mmseqs2_outputs')
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(mmseqs2_output_dir, exist_ok=True)
    
    # Download UniProt entries and extract capsid features
    print("=" * 80)
    print("Downloading UniProt entries and extracting capsid features")
    print("=" * 80)
    
    # Query for reviewed (SwissProt) Orthoretrovirinae Gag proteins
    # taxonomy_id:327045 = Orthoretrovirinae subfamily
    query = 'reviewed:true AND taxonomy_id:327045 AND gene:gag'
    json_dir = os.path.join(output_dir, 'orthoretrovirinae_gag_swissprot')
    
    total_downloaded = uniprot_utils.download_uniprot_jsons(
        query=query,
        batch_size=500,
        outdir=json_dir,
        verbose=False  # Use progress dots instead of per-file messages
    )
    
    all_hits = []
    json_files = [f for f in os.listdir(json_dir) if f.endswith('.json')]
    
    print(f"Processing {len(json_files)} JSON files...", end=" ", flush=True)
    
    files_with_hits = 0
    files_without_hits = 0
    
    for json_file in json_files:
        json_path = os.path.join(json_dir, json_file)
        
        with open(json_path, 'r') as f:
            entry = json.load(f)
        
        # Extract capsid features
        hits = uniprot_utils.extract_capsid_features_from_entry(entry)
        
        if hits:
            files_with_hits += 1
            all_hits.extend(hits)
        else:
            files_without_hits += 1
    
    print("done.")
    print(f"  Files with capsid features: {files_with_hits}")
    print(f"  Files without capsid features: {files_without_hits}")
    print(f"  Total capsid features found: {len(all_hits)}")
    
    # Save all hits metadata
    all_hits_file = os.path.join(output_dir, 'all_capsid_hits.json')
    metadata_utils.save_metadata_json(all_hits, all_hits_file, verbose=True)
    
    # Find unique sequences and save results
    print("\n" + "=" * 80)
    print("Finding unique sequences and saving results")
    print("=" * 80)
    
    unique_sequences = metadata_utils.aggregate_metadata_by_sequence(all_hits)
    n_unique = len(unique_sequences)
    print(f"Found {n_unique} unique sequences from {len(all_hits)} total hits")
    
    unique_sequences_file = os.path.join(output_dir, 'unique_capsid_sequences.json')
    unique_fasta = os.path.join(output_dir, 'unique_capsid_sequences.fasta')
    
    metadata_utils.save_metadata_json(unique_sequences, unique_sequences_file, verbose=False)
    fasta_utils.write_fasta_from_metadata_list(
        unique_sequences,
        unique_fasta,
        seq_key="sequence",
        id_key="primaryAccession"
    )
    print(f"Saved: {unique_sequences_file}, {unique_fasta}")
    
    # Cluster sequences and add cluster assignments
    print("\n" + "=" * 80)
    print("Clustering sequences with MMseqs2")
    print("=" * 80)
    
    cluster_output_prefix = os.path.join(mmseqs2_output_dir, "capsid_clusters")
    cluster_tmp_dir = os.path.join(output_dir, "mmseqs2_tmp")
    
    # Cluster sequences: 30% identity, 80% coverage, coverage of both query and target
    clustering_results = mmseqs2_utils.run_mmseqs_clustering(
        input_fasta=unique_fasta,
        output_prefix=cluster_output_prefix,
        tmp_dir=cluster_tmp_dir,
        min_seq_id=0.3,      # 30% minimum sequence identity
        coverage=0.8,         # 80% minimum coverage
        coverage_mode=0,      # Coverage of both query and target sequences
        remove_tmp_files=True
    )
    
    stats = clustering_results['stats']
    print(f"Clustering complete: {stats['n_clusters']} clusters from {stats['total_sequences']} sequences "
          f"(size range: {stats['min_cluster_size']}-{stats['max_cluster_size']}, "
          f"mean: {stats['mean_cluster_size']:.1f})")
    
    # Add cluster assignments and save updated sequences
    mmseqs2_utils.add_cluster_assignments(unique_sequences, clustering_results['output_files']['cluster_tsv'])
    metadata_utils.save_metadata_json(unique_sequences, unique_sequences_file, verbose=False)
    
    # Create TSV file
    tsv_file = os.path.join(output_dir, 'unique_capsid_sequences.tsv')
    metadata_utils.save_metadata_tsv(unique_sequences, tsv_file, verbose=True)
    
    # Summary
    print("\n" + "=" * 80)
    print("Summary")
    print("=" * 80)
    print(f"Total proteins processed: {len(json_files)}")
    print(f"Total capsid features found: {len(all_hits)}")
    print(f"Unique capsid sequences: {len(unique_sequences)}")
    print(f"Clusters: {stats['n_clusters']}")
    print(f"\nOutput files (in '{output_dir}' directory):")
    print(f"  - All hits: {all_hits_file}")
    print(f"  - Unique sequences JSON: {unique_sequences_file}")
    print(f"  - Unique sequences TSV: {tsv_file}")
    print(f"  - Unique sequences FASTA: {unique_fasta}")
    print(f"\nMMseqs2 output files (in '{mmseqs2_output_dir}' directory):")
    print(f"  - Cluster TSV: {clustering_results['output_files']['cluster_tsv']}")
    print(f"  - Representative sequences: {clustering_results['output_files']['rep_seq']}")
    print(f"  - All sequences per cluster: {clustering_results['output_files']['all_seqs']}")


if __name__ == "__main__":
    main()

# %%
