"""
UniProt API utilities for downloading and processing UniProt data.
"""
import os
import re
import json
import requests
from requests.adapters import HTTPAdapter, Retry


# Set up session with retries
_re_next_link = re.compile(r'<(.+)>; rel="next"')
_retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
_session = requests.Session()
_session.mount("https://", HTTPAdapter(max_retries=_retries))


def get_next_link(headers):
    """
    Extract the 'next' pagination link from HTTP Link header.
    
    UniProt API uses Link headers for pagination. This function extracts the URL
    for the next page of results from the 'rel="next"' link.
    
    Args:
        headers: Dictionary-like object containing HTTP response headers.
                Should contain a 'Link' header with pagination information.
    
    Returns:
        String URL for the next page if found, None otherwise.
    
    Example:
        Link header: '<https://rest.uniprot.org/uniprotkb/search?cursor=...>; rel="next"'
        Returns: 'https://rest.uniprot.org/uniprotkb/search?cursor=...'
    """
    if "Link" in headers:
        match = _re_next_link.match(headers["Link"])
        if match:
            return match.group(1)
    return None


def get_batch(batch_url, session=None):
    """
    Generator that yields batch responses from UniProt API with automatic pagination.
    
    This function handles pagination automatically by following 'next' links in the
    response headers. It continues fetching pages until all results are retrieved.
    
    Args:
        batch_url: Initial URL for the batch request. Should be a complete UniProt
                  REST API URL with query parameters. The function will follow
                  pagination links automatically.
        session: Optional requests.Session object. If not provided, uses the default
                session configured with retry logic. Useful for custom session
                configuration or testing.
    
    Yields:
        Tuple of (response, total) where:
            - response: requests.Response object containing the current page of results
            - total: Total number of results available (from X-Total-Results header),
                    or None if header is not present
    
    Raises:
        requests.HTTPError: If the HTTP request fails (raised by response.raise_for_status())
    
    Note:
        The generator will continue until no 'next' link is found in the response headers.
        Each yielded response contains a JSON payload with 'results' key containing
        the entries for that page.
    """
    if session is None:
        session = _session
    
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers.get("x-total-results", None)
        yield response, total
        batch_url = get_next_link(response.headers)


def download_uniprot_jsons(query='taxonomy_id:327046 AND gene:gag', batch_size=500, outdir='jsons', verbose=False):
    """
    Download UniProt entries matching a query as individual JSON files.
    
    This function queries the UniProt REST API, retrieves all matching entries,
    and saves each entry as a separate JSON file named by its primary accession.
    Handles pagination automatically to retrieve all results.
    
    Args:
        query: UniProt query string using UniProt query syntax.
              Examples:
              - 'taxonomy_id:327045 AND gene:gag' (Orthoretrovirinae Gag proteins)
              - 'reviewed:true AND organism:"Homo sapiens"'
              - 'accession:P12345'
              See UniProt documentation for full query syntax.
        batch_size: Number of results to retrieve per API request (default: 500).
                   Larger values reduce API calls but may hit rate limits.
                   Maximum recommended: 500.
        outdir: Output directory path for JSON files (default: 'jsons').
               Will be created if it doesn't exist.
               Each entry is saved as '{primaryAccession}.json'.
               If primaryAccession is missing, uses 'entry_{index}.json'.
        verbose: If True, prints each file saved. If False, shows progress dots (default: False).
    
    Returns:
        int: Total number of entries downloaded
    
    Note:
        Uses the default session with retry logic (5 retries with exponential backoff).
        Prints progress summary. If verbose=True, prints each file saved.
        If an entry lacks a primaryAccession, generates a default filename.
    """
    os.makedirs(outdir, exist_ok=True)
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": query,
        "format": "json",
        "size": batch_size
    }
    url = requests.Request('GET', base_url, params=params).prepare().url
    total_downloaded = 0
    for response, total in get_batch(url):
        data = response.json()
        for entry in data.get('results', []):
            accession = entry.get('primaryAccession', f'entry_{total_downloaded+1}')
            outpath = os.path.join(outdir, f"{accession}.json")
            with open(outpath, "w") as f:
                json.dump(entry, f)
            total_downloaded += 1
            if verbose:
                print(f"Saved {outpath}")
            elif total_downloaded % 50 == 0:
                print(".", end="", flush=True)
    if not verbose and total_downloaded > 0:
        print()  # New line after progress dots
    print(f"Downloaded {total_downloaded} entries as individual JSON files.")
    return total_downloaded


def extract_capsid_features_from_entry(entry):
    """
    Extract capsid protein features from a UniProt entry with complete metadata.
    
    This function searches through the features in a UniProt entry to find
    "Chain" features that are labeled as "capsid" (but not "nucleocapsid").
    Only complete features with exact start and end positions are extracted.
    
    Filtering criteria:
    - Feature type must be "Chain"
    - Description must contain "capsid" (case-insensitive)
    - Description must NOT contain "nucleocapsid"
    - Start and end positions must have modifier "EXACT" (complete features only)
    - End position must not extend to the end of the full protein sequence
    
    Args:
        entry: UniProt entry dictionary (typically from UniProt JSON API response).
              Should contain:
              - 'sequence': dict with 'value' key containing full protein sequence
              - 'primaryAccession': string
              - 'secondaryAccessions': list of strings
              - 'uniProtkbId': string
              - 'organism': dict with organism information
              - 'organismHosts': list of host dictionaries
              - 'features': list of feature dictionaries
    
    Returns:
        List of dictionaries, one for each capsid feature found. Each dictionary contains:
        - 'primaryAccession': Primary UniProt accession
        - 'secondaryAccessions': List of secondary accessions
        - 'uniProtkbId': UniProtKB identifier
        - 'organism': Dict with scientificName, commonName, taxonId, lineage
        - 'organismHosts': List of host dictionaries (scientificName, commonName, taxonId)
        - 'sequence': Extracted capsid sequence (substring of full protein sequence)
        - 'description': Feature description from UniProt
    
    Note:
        Returns empty list if no matching capsid features are found.
        Sequence positions are 1-based and inclusive (UniProt standard).
        The extracted sequence is the substring from start_pos-1 to end_pos (0-based indexing).
    """
    hits = []
    
    # Get full protein sequence (needed to extract capsid subsequence)
    full_sequence = entry.get("sequence", {}).get("value", "")
    
    # Extract basic entry metadata
    primary_accession = entry.get("primaryAccession", "")
    secondary_accessions = entry.get("secondaryAccessions", [])
    uniprotkb_id = entry.get("uniProtkbId", "")
    
    # Extract organism information (scientific name, common name, taxonomy)
    organism = entry.get("organism", {})
    organism_scientific_name = organism.get("scientificName", "")
    organism_common_name = organism.get("commonName", "")
    organism_taxon_id = organism.get("taxonId", "")
    organism_lineage = organism.get("lineage", [])
    
    # Extract organism hosts (viruses that infect this organism)
    organism_hosts = entry.get("organismHosts", [])
    hosts = []
    for host in organism_hosts:
        hosts.append({
            "scientificName": host.get("scientificName", ""),
            "commonName": host.get("commonName", ""),
            "taxonId": host.get("taxonId", "")
        })
    
    # Search through features to find capsid chains
    features = entry.get("features", [])
    for feature in features:
        # Only process "Chain" type features (not domains, regions, etc.)
        if feature.get("type") != "Chain":
            continue
        
        description = feature.get("description", "")
        desc_lower = description.lower()
        
        # Filter: must contain "capsid" but not "nucleocapsid" (we want capsid, not nucleocapsid)
        if "capsid" in desc_lower and "nucleocapsid" not in desc_lower:
            location = feature.get("location", {})
            start_info = location.get("start", {})
            end_info = location.get("end", {})
            
            # Only extract complete features (both start and end positions are EXACT)
            # This excludes partial or uncertain features
            start_modifier = start_info.get("modifier", "")
            end_modifier = end_info.get("modifier", "")
            
            if start_modifier == "EXACT" and end_modifier == "EXACT":
                start_pos = int(start_info.get("value", 0))
                end_pos = int(end_info.get("value", 0))
                
                # Skip if capsid extends to the end of the protein (likely incomplete annotation)
                if full_sequence and end_pos == len(full_sequence):
                    continue
                
                # Extract capsid sequence substring
                # UniProt uses 1-based indexing (inclusive), Python uses 0-based (exclusive end)
                # So convert: [start_pos, end_pos] -> [start_pos-1, end_pos)
                capsid_sequence = full_sequence[start_pos - 1:end_pos] if full_sequence else ""
                
                # Create hit record with all metadata
                hit = {
                    "primaryAccession": primary_accession,
                    "secondaryAccessions": secondary_accessions,
                    "uniProtkbId": uniprotkb_id,
                    "organism": {
                        "scientificName": organism_scientific_name,
                        "commonName": organism_common_name,
                        "taxonId": organism_taxon_id,
                        "lineage": organism_lineage
                    },
                    "organismHosts": hosts,
                    "sequence": capsid_sequence,  # Capsid sequence, not full protein
                    "description": description  # Description of the capsid feature
                }
                hits.append(hit)
    
    return hits
