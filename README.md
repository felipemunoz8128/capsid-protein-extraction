# Capsid Protein Extraction from Orthoretrovirinae Gag Proteins

A Python pipeline to extract and cluster capsid protein sequences from Orthoretrovirinae Gag proteins in UniProt SwissProt database.

## Overview

This project extracts capsid protein features from Orthoretrovirinae Gag proteins, identifies unique sequences, and clusters them using MMseqs2. The workflow:

1. Downloads UniProt SwissProt entries for Orthoretrovirinae Gag proteins
2. Extracts capsid features (complete, non-nucleocapsid) from each entry
3. Aggregates metadata by unique sequences
4. Clusters sequences using MMseqs2 (30% identity, 80% coverage)
5. Saves results in multiple formats (JSON, TSV, FASTA)

## Requirements

### Python Dependencies

- Python 3.7+
- `requests` - For UniProt API interactions
- `pandas` - For data manipulation

### External Tools

- **MMseqs2** - Must be installed and available in PATH as `mmseqs`

  Installation instructions:
  - Linux/Mac: `conda install -c conda-forge mmseqs2` or download from [MMseqs2 GitHub](https://github.com/soedinglab/MMseqs2)
  - See [MMseqs2 documentation](https://github.com/soedinglab/MMseqs2/wiki) for details

## Installation

1. Clone this repository:
```bash
git clone https://github.com/felipemunoz8128/capsid-protein-extraction.git
cd capsid-protein-extraction
```

2. Install Python dependencies:
```bash
pip install -r requirements.txt
```

3. Ensure MMseqs2 is installed and accessible:
```bash
mmseqs --version
```

## Usage

Run the main script:

```bash
python extract_capsid_proteins.py
```

The script will:
- Download UniProt entries matching the query: `reviewed:true AND taxonomy_id:327045 AND gene:gag`
- Extract capsid features from each entry
- Find unique sequences and aggregate metadata
- Cluster sequences using MMseqs2
- Save all results to the `outputs/` directory

### Output Files

All output files are saved in the `outputs/` directory:

- `all_capsid_hits.json` - All extracted capsid features with metadata
- `unique_capsid_sequences.json` - Unique sequences with aggregated metadata and cluster assignments
- `unique_capsid_sequences.tsv` - Tab-separated values file for easy viewing/analysis
- `unique_capsid_sequences.fasta` - FASTA format sequences
- `mmseqs2_outputs/` - MMseqs2 clustering results (cluster TSV, representative sequences, etc.)
- `orthoretrovirinae_gag_swissprot/` - Individual UniProt JSON entries downloaded from the API

**Note:** Example output files from a complete run are included in this repository in the `outputs/` directory. You can examine these files to see the expected output format before running the script yourself.

## Project Structure

```
capsid-protein-extraction/
├── extract_capsid_proteins.py  # Main workflow script
├── utils/                      # Utility modules
│   ├── __init__.py
│   ├── fasta_utils.py         # FASTA file operations
│   ├── metadata_utils.py      # Metadata aggregation and processing
│   ├── mmseqs2_utils.py      # MMseqs2 clustering operations
│   └── uniprot_utils.py       # UniProt API interactions
├── outputs/                    # Output files (example results included)
│   ├── all_capsid_hits.json
│   ├── unique_capsid_sequences.json
│   ├── unique_capsid_sequences.tsv
│   ├── unique_capsid_sequences.fasta
│   ├── mmseqs2_outputs/       # MMseqs2 clustering results
│   └── orthoretrovirinae_gag_swissprot/  # Downloaded UniProt entries
├── requirements.txt           # Python dependencies
└── README.md                  # This file
```

## Customization

You can modify the query in `extract_capsid_proteins.py` to extract different proteins:

```python
query = 'reviewed:true AND taxonomy_id:327045 AND gene:gag'
```

Clustering parameters can be adjusted in the `run_mmseqs_clustering()` call:

```python
clustering_results = mmseqs2_utils.run_mmseqs_clustering(
    input_fasta=unique_fasta,
    output_prefix=cluster_output_prefix,
    tmp_dir=cluster_tmp_dir,
    min_seq_id=0.3,      # Minimum sequence identity (0.0-1.0)
    coverage=0.8,         # Minimum coverage (0.0-1.0)
    coverage_mode=0,      # Coverage mode (0=both, 1=target, 2=query)
    remove_tmp_files=True
)
```


