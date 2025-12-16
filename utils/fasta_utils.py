"""
FASTA file utilities for reading, writing, and processing FASTA files.
"""


def _extract_label_from_uniprotkb_id(uniProtkb_id):
    """
    Extract organism label from uniProtkbId (part after underscore).
    
    Examples:
        "GAG_FIVWO" -> "FIVWO"
        ["GAG_FIVWO", "GAG_FIVCA"] -> "FIVWO,FIVCA"
        "GAG" -> "" (no underscore)
    
    Args:
        uniProtkb_id: String or list of strings containing UniProtKB IDs
    
    Returns:
        String with comma-separated labels, or empty string if no underscore found
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


def write_fasta_from_metadata_list(metadata_list, output_file, seq_key="sequence", id_key="primaryAccession", use_label=False):
    """
    Write a FASTA file from a list of metadata dictionaries.
    
    Each entry in the metadata list becomes a FASTA record with:
    - Header line: starts with '>' followed by the identifier
    - Sequence lines: the protein sequence (wrapped if needed by standard FASTA conventions)
    
    Args:
        metadata_list: List of dictionaries, each containing sequence and identifier info.
                      Each dictionary should contain at least the sequence and identifier keys.
        output_file: Output FASTA file path. Will be created or overwritten.
        seq_key: Key in each dict for the sequence (default: "sequence").
                If the value is a dict, looks for "value" key inside it.
                If the value is empty or missing, that entry is skipped.
        id_key: Key in each dict for the identifier (default: "primaryAccession").
                If the value is a list, all elements are joined with commas.
                If missing, generates a default identifier "seq_{index}".
        use_label: If True, uses the "label" field from metadata entries (if available),
                   or extracts label from uniProtkbId. The label field should already be
                   unique (duplicates handled as "label_primaryAccession" format).
                   Falls back to id_key if label is not available (default: False).
    
    Returns:
        None (writes directly to file)
    
    Note:
        Entries with empty sequences are skipped. The function handles both
        old format (sequence as dict with "value" key) and new format (sequence as direct string).
        When use_label=True, the function uses the "label" field from metadata entries,
        which should already be unique (duplicates are handled during metadata aggregation).
    """
    # Write FASTA file
    with open(output_file, 'w') as f:
        for i, entry in enumerate(metadata_list):
            # Get identifier
            if use_label:
                # Use label field if available (already handles duplicates)
                identifier = entry.get("label", "")
                # Fallback: extract label from uniProtkbId if label field not present
                if not identifier:
                    uniProtkb_id = entry.get("uniProtkbId", "")
                    identifier = _extract_label_from_uniprotkb_id(uniProtkb_id)
                # Final fallback to id_key if label extraction failed
                if not identifier:
                    identifier = entry.get(id_key, f"seq_{i}")
            else:
                identifier = entry.get(id_key, f"seq_{i}")
            
            if isinstance(identifier, list):
                # Join multiple accessions with commas
                identifier = ','.join(str(acc) for acc in identifier) if identifier else f"seq_{i}"
            
            # Get sequence (handle both old dict format and new direct format)
            seq_data = entry.get(seq_key, "")
            if isinstance(seq_data, dict):
                sequence = seq_data.get("value", "")
            else:
                sequence = str(seq_data) if seq_data else ""
            
            if sequence:
                f.write(f">{identifier}\n{sequence}\n")
