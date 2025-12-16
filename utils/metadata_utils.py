"""
Utilities for aggregating and processing metadata from UniProt entries.
"""
from collections import defaultdict


def _normalize_sequence(seq):
    """
    Normalize sequence format: handle both dict with "value" key and direct string.
    
    Args:
        seq: Sequence value (can be dict with "value" key or direct string)
    
    Returns:
        String sequence, or empty string if not found
    """
    if isinstance(seq, dict):
        return seq.get("value", "")
    return str(seq) if seq else ""


def aggregate_metadata_by_sequence(hits):
    """
    Group hits by sequence and aggregate metadata, then deduplicate by accession.
    
    This function performs two main operations:
    1. Groups hits by identical sequence and aggregates their metadata:
       - Fields with identical values across all hits are stored as single values
       - Fields with different values are stored as lists of unique values
       - Special handling for: primaryAccession, secondaryAccessions, uniProtkbId,
         organism fields, organismHosts, and description
    
    2. Deduplicates sequences that share the same primaryAccession(s):
       - If multiple sequences have the same primaryAccession(s), keeps only the longest sequence
       - This handles cases where the same protein has multiple sequence variants
    
    Args:
        hits: List of hit dictionaries, each containing:
              - "sequence": string or dict with "value" key containing the sequence
              - "primaryAccession": string or list of strings
              - "secondaryAccessions": list of strings
              - "uniProtkbId": string
              - "organism": dict with scientificName, commonName, taxonId, lineage
              - "organismHosts": list of host dicts
              - "description": string (optional)
              - Other metadata fields
    
    Returns:
        List of dictionaries with unique sequences and aggregated metadata.
        Each dictionary contains:
        - "sequence": the protein sequence (string)
        - "primaryAccession": string or list of strings (if multiple accessions)
        - "secondaryAccessions": list of unique secondary accessions
        - "uniProtkbId": string or list of strings (if multiple IDs)
        - "organism": dict with aggregated organism information
        - "organismHosts": list of unique host dictionaries
        - "description": string or list of strings (if multiple descriptions)
        - Other aggregated metadata fields
    """
    # Group hits by identical sequence (normalize format first)
    seq_to_hits = defaultdict(list)
    for hit in hits:
        seq = _normalize_sequence(hit["sequence"])
        seq_to_hits[seq].append(hit)
    
    # Aggregate metadata for each unique sequence
    unique_sequences_by_seq = []
    for sequence, group_hits in seq_to_hits.items():
        aggregated = {}
        
        # Aggregate metadata: if all hits have the same value, store as single value;
        # if values differ, store as list of unique values
        
        # Primary accession: single value if all hits share it, list if different
        accessions = [h["primaryAccession"] for h in group_hits]
        if len(set(accessions)) == 1:
            aggregated["primaryAccession"] = accessions[0]
        else:
            aggregated["primaryAccession"] = list(set(accessions))
        
        # Secondary accessions: collect all unique secondary accessions from all hits
        all_secondary = []
        for h in group_hits:
            all_secondary.extend(h.get("secondaryAccessions", []))
        aggregated["secondaryAccessions"] = list(set(all_secondary)) if all_secondary else []
        
        # UniProtKB ID
        uniprotkb_ids = [h["uniProtkbId"] for h in group_hits]
        if len(set(uniprotkb_ids)) == 1:
            aggregated["uniProtkbId"] = uniprotkb_ids[0]
        else:
            aggregated["uniProtkbId"] = list(set(uniprotkb_ids))
        
        # Organism - aggregate organism fields
        organisms_scientific = [h["organism"]["scientificName"] for h in group_hits]
        organisms_common = [h["organism"]["commonName"] for h in group_hits]
        organisms_taxon = [h["organism"]["taxonId"] for h in group_hits]
        
        aggregated["organism"] = {
            "scientificName": organisms_scientific[0] if len(set(organisms_scientific)) == 1 else list(set(organisms_scientific)),
            "commonName": organisms_common[0] if len(set(organisms_common)) == 1 else list(set(organisms_common)),
            "taxonId": organisms_taxon[0] if len(set(organisms_taxon)) == 1 else list(set(organisms_taxon)),
            "lineage": group_hits[0]["organism"]["lineage"]  # Lineage should be the same for same taxonomy
        }
        
        # Organism hosts: collect all unique hosts, deduplicating by taxonId
        all_hosts = []
        for h in group_hits:
            all_hosts.extend(h.get("organismHosts", []))
        # Remove duplicate hosts based on taxonId (same taxonId = same host)
        seen_taxon_ids = set()
        unique_hosts = []
        for host in all_hosts:
            taxon_id = host.get("taxonId")
            if taxon_id and taxon_id not in seen_taxon_ids:
                unique_hosts.append(host)
                seen_taxon_ids.add(taxon_id)
        aggregated["organismHosts"] = unique_hosts if unique_hosts else []
        
        # Description
        descriptions = [h.get("description", "") for h in group_hits if "description" in h]
        if descriptions:
            if len(set(descriptions)) == 1:
                aggregated["description"] = descriptions[0]
            else:
                aggregated["description"] = list(set(descriptions))
        
        # Sequence
        aggregated["sequence"] = sequence
        
        unique_sequences_by_seq.append(aggregated)
    
    # Second pass: deduplicate sequences that share the same primaryAccession(s)
    # If multiple sequences have the same accession(s), keep only the longest sequence
    # This handles cases where the same protein has multiple sequence variants
    accession_to_entry = {}  # normalized_accession_tuple -> (entry, sequence_length)
    
    for entry in unique_sequences_by_seq:
        # Normalize primaryAccession to a tuple for comparison (handles both single and list formats)
        acc = entry.get("primaryAccession", "")
        if isinstance(acc, list):
            # Sort and deduplicate, filter out empty strings for consistent comparison
            acc_list = sorted(set(str(a) for a in acc if a))
            acc_tuple = tuple(acc_list) if acc_list else ("",)
        else:
            acc_tuple = (str(acc),) if acc else ("",)
        
        # Extract sequence length (normalize format first)
        seq = _normalize_sequence(entry["sequence"])
        seq_length = len(seq)
        
        # Keep the longest sequence for each unique set of accessions
        if acc_tuple not in accession_to_entry:
            accession_to_entry[acc_tuple] = (entry, seq_length)
        else:
            existing_entry, existing_length = accession_to_entry[acc_tuple]
            if seq_length > existing_length:
                # Replace with longer sequence
                accession_to_entry[acc_tuple] = (entry, seq_length)
            elif seq_length == existing_length:
                # Same length - keep the first one encountered (arbitrary but deterministic)
                pass
    
    # Return only the entries we kept
    unique_sequences = [entry for entry, _ in accession_to_entry.values()]
    
    return unique_sequences


def save_metadata_json(metadata_list, output_file, verbose=True):
    """
    Save metadata list to a JSON file with pretty-printed formatting.
    
    The output JSON file is formatted with 2-space indentation for readability.
    Optionally prints a confirmation message with the number of entries saved.
    
    Args:
        metadata_list: List of metadata dictionaries to save. Each dictionary
                      will be serialized as a JSON object in the output array.
        output_file: Output file path. Will be created or overwritten.
                     Should have .json extension for clarity.
        verbose: If True, prints confirmation message (default: True).
    
    Returns:
        None (writes directly to file)
    
    Note:
        Uses json.dump() with indent=2 for human-readable formatting.
        All entries are saved in a single JSON array.
    """
    import json
    
    with open(output_file, "w") as f:
        json.dump(metadata_list, f, indent=2)
    if verbose:
        print(f"Saved {len(metadata_list)} entries to {output_file}")


def save_metadata_tsv(metadata_list, output_file, verbose=True):
    """
    Save metadata list to a TSV (tab-separated values) file with one row per entry.
    
    Creates a tabular representation of the metadata with the following columns:
    - label: Extracted from uniProtkbId (part after underscore, e.g., "FIVWO" from "GAG_FIVWO")
    - primaryAccession: Primary UniProt accession(s), comma-separated if multiple
    - secondaryAccessions: Secondary accessions, comma-separated
    - uniProtkbId: UniProtKB identifier(s), comma-separated if multiple
    - scientificName: Organism scientific name(s), comma-separated if multiple
    - commonName: Organism common name(s), comma-separated if multiple
    - taxonId: Organism taxonomy ID(s), comma-separated if multiple
    - genus: Extracted from organism lineage (value following "Orthoretrovirinae")
    - sequence: Protein sequence string
    - description: Feature description(s), comma-separated if multiple
    - cluster_id: MMseqs2 cluster assignment (if available)
    
    Args:
        metadata_list: List of metadata dictionaries. Each dictionary becomes one row.
                      Missing fields are represented as empty strings.
        output_file: Output TSV file path. Will be created or overwritten.
                    Should have .tsv extension for clarity.
        verbose: If True, prints confirmation message (default: True).
    
    Returns:
        None (writes directly to file)
    
    Note:
        - List values are joined with commas
        - Empty/None values are represented as empty strings
        - The file uses UTF-8 encoding
        - First row contains column headers
    """
    import csv
    
    if not metadata_list:
        if verbose:
            print("No data to save")
        return
    
    def extract_organism_labels(uniProtkb_id):
        """
        Extract organism labels from uniProtkbId by taking the part after the underscore.
        
        Examples:
            "GAG_FIVWO" -> "FIVWO"
            ["GAG_FIVWO", "GAG_FIVCA"] -> "FIVWO,FIVCA"
            "GAG" -> "" (no underscore)
        
        Args:
            uniProtkb_id: String or list of strings containing UniProtKB IDs
        
        Returns:
            String with comma-separated labels (part after underscore), or empty string if no underscore found
        """
        if isinstance(uniProtkb_id, list):
            labels = []
            for uid in uniProtkb_id:
                if '_' in str(uid):
                    labels.append(str(uid).split('_', 1)[1])
            return ','.join(labels) if labels else ''
        elif isinstance(uniProtkb_id, str) and '_' in uniProtkb_id:
            return uniProtkb_id.split('_', 1)[1]
        return ''
    
    def get_genus(entry):
        """
        Extract genus name from organism lineage.
        
        Looks for "Orthoretrovirinae" in the lineage list and returns the
        value immediately following it, which should be the genus name.
        
        Args:
            entry: Metadata dictionary containing 'organism' key with 'lineage' list
        
        Returns:
            String genus name if found, empty string otherwise
        
        Example:
            lineage: ["Viruses", "Riboviria", ..., "Orthoretrovirinae", "Lentivirus", ...]
            Returns: "Lentivirus"
        """
        lineage = entry.get('organism', {}).get('lineage', [])
        if isinstance(lineage, list):
            try:
                idx = lineage.index('Orthoretrovirinae')
                if idx + 1 < len(lineage):
                    return lineage[idx + 1]
            except ValueError:
                pass
        return ''
    
    def format_value(value):
        """
        Format a value for TSV output by converting to string and handling lists.
        
        Args:
            value: Value to format (can be string, list, None, or other types)
        
        Returns:
            String representation:
            - Lists: comma-separated string of non-empty values
            - None: empty string
            - Other types: string conversion
        """
        if isinstance(value, list):
            return ','.join(str(v) for v in value if v)
        elif value is None:
            return ''
        else:
            return str(value)
    
    # Define column order
    columns = ['label', 'primaryAccession', 'secondaryAccessions', 'uniProtkbId',
               'scientificName', 'commonName', 'taxonId', 'genus', 
               'sequence', 'description', 'cluster_id']
    
    # Write TSV file
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        
        # Write header
        writer.writerow(columns)
        
        # Write data rows
        for entry in metadata_list:
            row = []
            
            # Extract organism labels from uniProtkbId(s)
            uniProtkb_id = entry.get('uniProtkbId', '')
            label = extract_organism_labels(uniProtkb_id)
            row.append(label)
            
            # Primary accession
            row.append(format_value(entry.get('primaryAccession', '')))
            
            # Secondary accessions
            row.append(format_value(entry.get('secondaryAccessions', [])))
            
            # UniProtKB ID
            row.append(format_value(entry.get('uniProtkbId', '')))
            
            # Organism fields
            organism = entry.get('organism', {})
            row.append(format_value(organism.get('scientificName', '')))
            row.append(format_value(organism.get('commonName', '')))
            row.append(format_value(organism.get('taxonId', '')))
            
            # Genus
            row.append(get_genus(entry))
            
            # Sequence
            seq = _normalize_sequence(entry.get('sequence', ''))
            row.append(str(seq) if seq else '')
            
            # Description
            row.append(format_value(entry.get('description', '')))
            
            # Cluster ID
            row.append(format_value(entry.get('cluster_id', '')))
            
            writer.writerow(row)
    
    if verbose:
        print(f"Saved {len(metadata_list)} entries to {output_file}")

