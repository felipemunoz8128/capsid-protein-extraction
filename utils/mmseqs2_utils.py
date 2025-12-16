"""
Utilities for running MMseqs2 clustering and parsing results.
"""
import os
import shutil
import subprocess
import pandas as pd
from utils.fasta_utils import _extract_label_from_uniprotkb_id


def run_mmseqs_clustering(
    input_fasta: str,
    output_prefix: str,
    tmp_dir: str,
    min_seq_id: float = 0.3,
    coverage: float = 0.8,
    coverage_mode: int = 0,
    remove_tmp_files: bool = True
) -> dict:
    """
    Run MMseqs2 easy-cluster on a protein dataset to group similar sequences.
    
    This function executes MMseqs2's easy-cluster command, which performs sequence
    clustering based on sequence identity and coverage thresholds. It automatically
    handles pagination and temporary file management.
    
    Args:
        input_fasta: Path to input FASTA file containing protein sequences to cluster.
                    Each sequence header should uniquely identify the sequence.
        output_prefix: Prefix for output files. The following files will be created:
                      - {output_prefix}_cluster.tsv: Cluster assignments (rep_id, member_id)
                      - {output_prefix}_rep_seq.fasta: Representative sequences for each cluster
                      - {output_prefix}_all_seqs.fasta: All sequences grouped by cluster
        tmp_dir: Directory for temporary MMseqs2 files. Will be created if it doesn't exist.
                Can be safely deleted after clustering completes if remove_tmp_files=True.
        min_seq_id: Minimum sequence identity threshold (0.0-1.0). Sequences must share
                   at least this fraction of identical residues to be clustered together.
                   Default: 0.3 (30% identity).
        coverage: Minimum coverage threshold (0.0-1.0). The alignment must cover at least
                 this fraction of the sequence length. Default: 0.8 (80% coverage).
        coverage_mode: Coverage calculation mode:
                      - 0: Coverage of both query and target sequences (default)
                      - 1: Coverage of target sequence only
                      - 2: Coverage of query sequence only
        remove_tmp_files: If True, removes temporary MMseqs2 files after clustering
                         completes. If False, keeps them for debugging. Default: True.
    
    Returns:
        Dictionary with the following keys:
            - 'output_files': dict with keys 'cluster_tsv', 'rep_seq', 'all_seqs'
                             containing paths to the respective output files
            - 'stats': dict with clustering statistics:
                      - 'n_clusters': number of clusters found
                      - 'total_sequences': total number of sequences clustered
                      - 'min_cluster_size': size of smallest cluster
                      - 'max_cluster_size': size of largest cluster
                      - 'mean_cluster_size': average cluster size
                      - 'median_cluster_size': median cluster size
            - 'cluster_sizes': list of cluster sizes (one value per cluster)
    
    Raises:
        FileNotFoundError: If expected output files are not created after clustering
        subprocess.CalledProcessError: If MMseqs2 command fails
        Exception: For other unexpected errors during clustering
    
    Note:
        Requires MMseqs2 to be installed and available in PATH as 'mmseqs'.
        The function prints progress messages and error details if clustering fails.
    """
    # Create temporary directory if it doesn't exist
    os.makedirs(tmp_dir, exist_ok=True)
    
    # Construct the MMseqs2 command
    cmd = [
        "mmseqs",
        "easy-cluster",
        input_fasta,
        output_prefix,
        tmp_dir,
        "--min-seq-id", str(min_seq_id),
        "-c", str(coverage),
        "--cov-mode", str(coverage_mode)
    ]
    
    if remove_tmp_files:
        cmd.append("--remove-tmp-files")
    
    try:
        # Run the command
        print(f"Running MMseqs2 clustering with min-seq-id={min_seq_id}, coverage={coverage}, cov-mode={coverage_mode}...")
        subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        
        # Get output file paths
        output_files = {
            'cluster_tsv': f"{output_prefix}_cluster.tsv",
            'rep_seq': f"{output_prefix}_rep_seq.fasta",
            'all_seqs': f"{output_prefix}_all_seqs.fasta"
        }
        
        # Verify output files exist
        for file_type, file_path in output_files.items():
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Expected output file {file_path} not found")
        
        # Parse cluster assignments and calculate statistics
        n_clusters, cluster_sizes, stats = parse_cluster_tsv(output_files['cluster_tsv'])
        
        # Explicitly remove temporary directory if requested
        # (MMseqs2 --remove-tmp-files may not remove the directory itself)
        if remove_tmp_files and os.path.exists(tmp_dir):
            try:
                shutil.rmtree(tmp_dir)
            except Exception as e:
                # Don't fail if cleanup fails, just warn
                print(f"Warning: Could not remove temporary directory {tmp_dir}: {e}")
        
        return {
            'output_files': output_files,
            'stats': stats,
            'cluster_sizes': cluster_sizes
        }
        
    except subprocess.CalledProcessError as e:
        print(f"Error running MMseqs2: {e}")
        if e.stderr:
            print(f"Error output: {e.stderr[:500]}")
        raise
    except Exception as e:
        print(f"Unexpected error: {e}")
        raise


def _parse_cluster_members(cluster_tsv: str):
    """
    Parse MMseqs2 cluster TSV file and return cluster membership dictionary.
    
    Helper function that extracts the core parsing logic shared by parse_cluster_tsv
    and get_cluster_dataframe.
    
    Args:
        cluster_tsv: Path to cluster TSV file created by MMseqs2 easy-cluster command.
    
    Returns:
        Dictionary mapping representative_id -> list of member_ids (including representative)
    """
    cluster_members = {}
    
    with open(cluster_tsv, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                rep_id = parts[0]      # Representative sequence ID for this cluster
                member_id = parts[1]    # Member sequence ID (same as rep_id for singletons)
                
                if rep_id not in cluster_members:
                    cluster_members[rep_id] = []
                
                # Add member if not already present (avoid duplicate entries)
                if member_id not in cluster_members[rep_id]:
                    cluster_members[rep_id].append(member_id)
    
    return cluster_members


def parse_cluster_tsv(cluster_tsv: str):
    """
    Parse MMseqs2 cluster TSV file and extract cluster statistics.
    
    The TSV format is: representative_id<TAB>member_id
    Each line represents a sequence in a cluster, where:
    - First column: representative sequence ID for that cluster
    - Second column: member sequence ID (can be the same as representative for singletons)
    
    The same representative ID can appear on multiple lines if the cluster has multiple members.
    
    Args:
        cluster_tsv: Path to cluster TSV file created by MMseqs2 easy-cluster command.
                    File should be tab-separated with two columns per line.
    
    Returns:
        Tuple containing:
            - n_clusters (int): Number of unique clusters (number of unique representatives)
            - cluster_sizes (list): List of cluster sizes, one integer per cluster
            - stats (dict): Dictionary with clustering statistics:
                           - 'n_clusters': number of clusters
                           - 'total_sequences': total number of sequences
                           - 'min_cluster_size': minimum cluster size
                           - 'max_cluster_size': maximum cluster size
                           - 'mean_cluster_size': average cluster size
                           - 'median_cluster_size': median cluster size
    
    Note:
        Empty lines are skipped. Duplicate entries (same rep_id, member_id pair)
        are ignored (only counted once per cluster).
    """
    cluster_members = _parse_cluster_members(cluster_tsv)
    
    n_clusters = len(cluster_members)
    cluster_sizes = [len(members) for members in cluster_members.values()]
    
    stats = {
        "n_clusters": n_clusters,
        "total_sequences": sum(cluster_sizes),
        "min_cluster_size": min(cluster_sizes) if cluster_sizes else 0,
        "max_cluster_size": max(cluster_sizes) if cluster_sizes else 0,
        "mean_cluster_size": sum(cluster_sizes) / len(cluster_sizes) if cluster_sizes else 0,
        "median_cluster_size": sorted(cluster_sizes)[len(cluster_sizes) // 2] if cluster_sizes else 0
    }
    
    return n_clusters, cluster_sizes, stats


def get_cluster_dataframe(cluster_tsv: str) -> pd.DataFrame:
    """
    Parse MMseqs2 cluster TSV file into a pandas DataFrame for easier manipulation.
    
    Creates a DataFrame where each row represents a sequence-cluster assignment.
    The DataFrame includes cluster size information for each row.
    
    Args:
        cluster_tsv: Path to cluster TSV file created by MMseqs2 easy-cluster command.
                    File should be tab-separated with two columns per line.
    
    Returns:
        pandas.DataFrame with columns:
            - 'representative_id': Representative sequence ID for the cluster
            - 'member_id': Member sequence ID (can be same as representative)
            - 'cluster_size': Number of sequences in this cluster (integer)
        
        Each row represents one sequence's cluster assignment. If a cluster has
        N members, there will be N rows with the same representative_id and cluster_size.
    
    Note:
        Empty lines are skipped. Duplicate entries are ignored.
        Requires pandas to be installed.
        Internally reuses shared parsing logic to avoid code duplication.
    """
    # Reuse shared parsing logic
    cluster_members = _parse_cluster_members(cluster_tsv)
    
    # Build DataFrame
    rows = []
    for rep_id, members in cluster_members.items():
        for member_id in members:
            rows.append({
                'representative_id': rep_id,
                'member_id': member_id,
                'cluster_size': len(members)
            })
    
    return pd.DataFrame(rows)


def add_cluster_assignments(metadata_list, cluster_tsv):
    """
    Add cluster_id field to metadata entries based on MMseqs2 cluster assignments.
    
    This function matches sequences in the metadata list with their cluster assignments
    from the MMseqs2 output. Cluster IDs are assigned sequentially starting from 1,
    based on the order representatives appear in the cluster TSV file.
    
    The function modifies the metadata_list in place by adding a 'cluster_id' field
    to each entry. Entries that cannot be matched receive cluster_id=None.
    
    Args:
        metadata_list: List of metadata dictionaries to update. Each dictionary should
                     have a 'label' field (added during metadata aggregation) that matches
                     the FASTA headers used for clustering. Modified in place.
        cluster_tsv: Path to MMseqs2 cluster TSV file. Should be the same file used
                    to cluster the sequences (typically ends with '_cluster.tsv').
    
    Returns:
        Dictionary mapping sequence identifiers (as they appear in the cluster TSV)
        to cluster IDs (integers starting from 1). This can be used for custom
        lookup logic if needed.
    
    Note:
        The function performs multiple lookup strategies to match sequences:
        1. Uses the 'label' field from metadata entries (preferred, already unique)
        2. Falls back to extracting label from uniProtkbId if label field not present
        3. Falls back to primaryAccession if label extraction fails
        4. Tries sorted comma-separated format if lookup_key contains commas
        5. Tries individual labels/accessions as last resort
        If no match is found, cluster_id is set to None.
        
        Cluster IDs are assigned based on representative sequence order, not cluster size.
        The function assumes the FASTA file was written with labels (use_label=True),
        and that the 'label' field in metadata entries matches the FASTA headers.
    """
    # Load cluster assignments
    cluster_df = get_cluster_dataframe(cluster_tsv)
    
    # Build mapping from sequence identifier (as it appears in FASTA) to cluster_id
    # Cluster IDs are assigned sequentially (1, 2, 3, ...) based on representative order
    sequence_to_cluster_id = {}
    representative_to_cluster_num = {}
    
    for idx, row in cluster_df.iterrows():
        rep_id = row['representative_id']
        member_id = row['member_id']
        
        # Assign sequential cluster number to each unique representative (first time we see it)
        if rep_id not in representative_to_cluster_num:
            representative_to_cluster_num[rep_id] = len(representative_to_cluster_num) + 1
        
        # All members of a cluster get the same cluster_id as their representative
        cluster_id = representative_to_cluster_num[rep_id]
        sequence_to_cluster_id[member_id] = cluster_id
    
    # Add cluster_id to each metadata entry by matching sequence identifiers
    # FASTA file uses the "label" field from metadata entries, which is already unique
    for entry in metadata_list:
        # Use label field if available (already handles duplicates)
        lookup_key = entry.get('label', '')
        
        # Fallback: extract label from uniProtkbId if label field not present
        if not lookup_key:
            uniProtkb_id = entry.get('uniProtkbId', '')
            lookup_key = _extract_label_from_uniprotkb_id(uniProtkb_id)
        
        # Final fallback to primaryAccession if label extraction failed
        if not lookup_key:
            acc = entry.get('primaryAccession', '')
            if isinstance(acc, list):
                lookup_key = ','.join(str(a) for a in acc if a)
            else:
                lookup_key = str(acc) if acc else ""
        
        # Try exact match
        cluster_id = sequence_to_cluster_id.get(lookup_key)
        
        if cluster_id is not None:
            entry['cluster_id'] = cluster_id
        else:
            # Strategy 2: If not found and lookup_key contains commas, try sorted version
            if ',' in lookup_key:
                sorted_key = ','.join(sorted(lookup_key.split(',')))
                cluster_id = sequence_to_cluster_id.get(sorted_key)
                if cluster_id is not None:
                    entry['cluster_id'] = cluster_id
                else:
                    # Strategy 3: Try individual labels/accessions (last resort)
                    parts = lookup_key.split(',')
                    if len(parts) > 0:
                        # Try each part individually
                        found = False
                        for part in parts:
                            cluster_id = sequence_to_cluster_id.get(part.strip())
                            if cluster_id is not None:
                                entry['cluster_id'] = cluster_id
                                found = True
                                break
                        if not found:
                            entry['cluster_id'] = None
                    else:
                        entry['cluster_id'] = None
            else:
                entry['cluster_id'] = None
    
    return sequence_to_cluster_id

