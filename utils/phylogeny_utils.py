"""
Utilities for phylogenetic analysis preparation: multiple sequence alignment and trimming.
"""
import os
import shutil
import subprocess
from typing import Optional


def run_mafft_alignment(
    input_fasta: str,
    output_fasta: str,
    algorithm: str = "linsi",
    maxiterate: Optional[int] = None,
    op: Optional[float] = None,
    ep: Optional[float] = None,
    threads: Optional[int] = None,
    log_file: Optional[str] = None,
    verbose: bool = True
) -> dict:
    """
    Run MAFFT multiple sequence alignment on a FASTA file.
    
    Args:
        input_fasta: Path to input FASTA file containing sequences to align
        output_fasta: Path to output aligned FASTA file
        algorithm: Algorithm to use:
                  - "auto": Automatically select (default for MAFFT)
                  - "linsi": High accuracy, local pairwise (recommended for proteins)
                  - "einsi": High accuracy, genafpair
                  - "ginsi": High accuracy, globalpair
                  - "fast": High speed (--retree 1)
        maxiterate: Maximum number of iterative refinement (overrides algorithm default)
        op: Gap opening penalty (default: 1.53)
        ep: Offset/gap extension penalty (default: 0.0)
        threads: Number of threads (-1 for auto, None to use all available CPUs)
        log_file: Path to log file for saving full MAFFT output (default: None, no log saved)
        verbose: Print progress messages (default: True)
    
    Returns:
        Dictionary with:
        - 'success' (bool): Whether alignment succeeded
        - 'output_file' (str): Path to output file
        - 'log_file' (str, optional): Path to log file if created
        - 'error' (str, optional): Error message if failed
    """
    # Validate input file exists
    if not os.path.exists(input_fasta):
        raise FileNotFoundError(f"Input FASTA file not found: {input_fasta}")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_fasta)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Build command
    cmd = ["mafft"]
    
    # Algorithm selection
    if algorithm == "auto":
        cmd.append("--auto")
    elif algorithm == "fast":
        cmd.extend(["--retree", "1"])
    elif algorithm == "linsi":
        cmd.extend(["--maxiterate", "1000", "--localpair"])
    elif algorithm == "einsi":
        cmd.extend(["--maxiterate", "1000", "--genafpair"])
    elif algorithm == "ginsi":
        cmd.extend(["--maxiterate", "1000", "--globalpair"])
    else:
        # Use algorithm as-is (user can pass custom flags)
        if algorithm:
            cmd.append(algorithm)
    
    # Override maxiterate if specified
    if maxiterate is not None:
        # Remove any existing --maxiterate
        cmd = [c for c in cmd if not c.startswith("--maxiterate")]
        cmd.extend(["--maxiterate", str(maxiterate)])
    
    # Gap penalties
    if op is not None:
        cmd.extend(["--op", str(op)])
    if ep is not None:
        cmd.extend(["--ep", str(ep)])
    
    # Threads (default to -1 for auto if not specified)
    if threads is None:
        threads = -1
    cmd.extend(["--thread", str(threads)])
    
    # Input file
    cmd.append(input_fasta)
    
    if verbose:
        print(f"Running MAFFT alignment (algorithm: {algorithm})...")
    
    # Run command
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Write alignment output to file
        with open(output_fasta, 'w') as f:
            f.write(result.stdout)
        
        # Save log file if requested
        if log_file:
            with open(log_file, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("MAFFT Alignment Log\n")
                f.write("=" * 80 + "\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(f"Algorithm: {algorithm}\n")
                f.write(f"Input: {input_fasta}\n")
                f.write(f"Output: {output_fasta}\n")
                f.write("\n" + "=" * 80 + "\n")
                f.write("STDOUT:\n")
                f.write("=" * 80 + "\n")
                if result.stderr:
                    f.write(result.stderr)
                f.write("\n" + "=" * 80 + "\n")
                f.write("STDERR:\n")
                f.write("=" * 80 + "\n")
                if result.stderr:
                    f.write(result.stderr)
        
        if verbose:
            print(f"Alignment completed successfully: {output_fasta}")
            if log_file:
                print(f"Log saved: {log_file}")
        
        return_dict = {
            "success": True,
            "output_file": output_fasta
        }
        if log_file:
            return_dict["log_file"] = log_file
        
        return return_dict
        
    except subprocess.CalledProcessError as e:
        error_msg = f"MAFFT alignment failed with return code {e.returncode}"
        
        # Save error to log file if requested
        if log_file:
            with open(log_file, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("MAFFT Alignment Log - ERROR\n")
                f.write("=" * 80 + "\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(f"Return code: {e.returncode}\n")
                f.write(f"STDOUT:\n{e.stdout}\n")
                f.write(f"STDERR:\n{e.stderr}\n")
        
        if verbose:
            print(f"Error: {error_msg}")
            if e.stderr:
                print(f"Error details: {e.stderr[:500]}")
            if log_file:
                print(f"Error log saved: {log_file}")
        
        return {
            "success": False,
            "error": error_msg
        }
    except FileNotFoundError:
        error_msg = "MAFFT not found. Please ensure MAFFT is installed and in your PATH."
        if verbose:
            print(f"Error: {error_msg}")
        return {
            "success": False,
            "error": error_msg
        }


def run_clipkit_trimming(
    input_alignment: str,
    output_alignment: str,
    mode: str = "smart-gap",
    gaps: float = 0.9,
    sequence_type: Optional[str] = None,
    threads: int = 1,
    clipkit_path: Optional[str] = None,
    log_file: Optional[str] = None,
    verbose: bool = True
) -> dict:
    """
    Run ClipKIT to trim a multiple sequence alignment.
    
    Args:
        input_alignment: Path to input alignment file (MAFFT output)
        output_alignment: Path to output trimmed alignment file
        mode: Trimming mode:
              - "smart-gap": Dynamic determination of gaps threshold (default, recommended)
              - "gappy": Trim sites > gaps threshold
              - "kpic": Keep parsimony informative and constant sites
              - "kpic-smart-gap": Combination of kpic and smart-gap
              - "kpi": Keep only parsimony informative sites
        gaps: Gaps threshold (0-1, default: 0.9). Ignored for smart-gap mode.
        sequence_type: Sequence type: "aa" (amino acid) or "nt" (nucleotide).
                      If None, auto-detects from alignment.
        threads: Number of threads for parallel processing (default: 1)
        clipkit_path: Path to clipkit executable. If None, attempts to find automatically.
        log_file: Path to log file for saving full ClipKIT output (default: None, no log saved)
        verbose: Print progress messages (default: True)
    
    Returns:
        Dictionary with:
        - 'success' (bool): Whether trimming succeeded
        - 'output_file' (str): Path to output file
        - 'log_file' (str, optional): Path to log file if created
        - 'error' (str, optional): Error message if failed
    """
    # Validate input file exists
    if not os.path.exists(input_alignment):
        raise FileNotFoundError(f"Input alignment file not found: {input_alignment}")
    
    # Find clipkit executable in PATH
    if clipkit_path is None:
        clipkit_path = shutil.which("clipkit")
    
    if clipkit_path is None:
        error_msg = "ClipKIT not found in PATH. Please ensure the conda environment is activated: conda activate capsid-protein-extraction"
        if verbose:
            print(f"Error: {error_msg}")
        return {
            "success": False,
            "error": error_msg
        }
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_alignment)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Build command
    cmd = [clipkit_path, input_alignment, "-o", output_alignment]
    
    # Mode
    cmd.extend(["-m", mode])
    
    # Gaps threshold
    cmd.extend(["-g", str(gaps)])
    
    # Sequence type
    if sequence_type:
        cmd.extend(["-s", sequence_type])
    
    # Threads
    if threads > 1:
        cmd.extend(["-t", str(threads)])
    
    # Don't use -q flag so we can capture full output for log file
    # We'll suppress console output by not printing stdout
    
    if verbose:
        print(f"Running ClipKIT trimming (mode: {mode})...")
    
    # Run command
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=True
        )
        
        # Save log file if requested (contains full ClipKIT output)
        if log_file:
            with open(log_file, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("ClipKIT Trimming Log\n")
                f.write("=" * 80 + "\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(f"Mode: {mode}\n")
                f.write(f"Gaps threshold: {gaps}\n")
                f.write(f"Input: {input_alignment}\n")
                f.write(f"Output: {output_alignment}\n")
                f.write("\n" + "=" * 80 + "\n")
                f.write("STDOUT:\n")
                f.write("=" * 80 + "\n")
                if result.stdout:
                    f.write(result.stdout)
                f.write("\n" + "=" * 80 + "\n")
                f.write("STDERR:\n")
                f.write("=" * 80 + "\n")
                if result.stderr:
                    f.write(result.stderr)
        
        if verbose:
            print(f"Trimming completed successfully: {output_alignment}")
            if log_file:
                print(f"Log saved: {log_file}")
        
        return_dict = {
            "success": True,
            "output_file": output_alignment
        }
        if log_file:
            return_dict["log_file"] = log_file
        
        return return_dict
        
    except subprocess.CalledProcessError as e:
        error_msg = f"ClipKIT trimming failed with return code {e.returncode}"
        
        # Save error to log file if requested
        if log_file:
            with open(log_file, 'w') as f:
                f.write("=" * 80 + "\n")
                f.write("ClipKIT Trimming Log - ERROR\n")
                f.write("=" * 80 + "\n")
                f.write(f"Command: {' '.join(cmd)}\n")
                f.write(f"Return code: {e.returncode}\n")
                f.write(f"STDOUT:\n{e.stdout}\n")
                f.write(f"STDERR:\n{e.stderr}\n")
        
        if verbose:
            print(f"Error: {error_msg}")
            if e.stderr:
                print(f"Error details: {e.stderr[:500]}")
            if log_file:
                print(f"Error log saved: {log_file}")
        
        return {
            "success": False,
            "error": error_msg
        }
    except FileNotFoundError:
        error_msg = f"ClipKIT executable not found at: {clipkit_path}"
        if verbose:
            print(f"Error: {error_msg}")
        return {
            "success": False,
            "error": error_msg
        }


def prepare_sequences_for_phylogeny(
    input_fasta: str,
    output_dir: str,
    mafft_algorithm: str = "linsi",
    clipkit_mode: str = "smart-gap",
    save_logs: bool = True,
    verbose: bool = True
) -> dict:
    """
    Complete workflow to prepare sequences for phylogenetic tree building.
    
    This function runs:
    1. MAFFT multiple sequence alignment
    2. ClipKIT alignment trimming
    
    The output trimmed alignment is ready for tree inference with IQTREE3.
    
    Args:
        input_fasta: Path to input FASTA file with sequences to align
        output_dir: Directory for output files
        mafft_algorithm: MAFFT algorithm to use (default: "linsi" for high accuracy)
        clipkit_mode: ClipKIT trimming mode (default: "smart-gap")
        save_logs: If True, save full command output to log files (default: True)
        verbose: Print progress messages (default: True)
    
    Returns:
        Dictionary with:
        - 'success' (bool): Whether all steps succeeded
        - 'aligned_file' (str): Path to MAFFT-aligned sequences
        - 'trimmed_file' (str): Path to ClipKIT-trimmed alignment (ready for IQTREE3)
        - 'mafft_log' (str, optional): Path to MAFFT log file if saved
        - 'clipkit_log' (str, optional): Path to ClipKIT log file if saved
        - 'error' (str, optional): Error message if any step failed
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # File paths
    aligned_file = os.path.join(output_dir, "aligned_sequences.fasta")
    trimmed_file = os.path.join(output_dir, "aligned_sequences_trimmed.fasta")
    mafft_log_file = os.path.join(output_dir, "mafft_log.txt") if save_logs else None
    clipkit_log_file = os.path.join(output_dir, "clipkit_log.txt") if save_logs else None
    
    # Step 1: MAFFT alignment
    if verbose:
        print(f"Aligning {os.path.basename(input_fasta)} with MAFFT...")
    
    align_result = run_mafft_alignment(
        input_fasta=input_fasta,
        output_fasta=aligned_file,
        algorithm=mafft_algorithm,
        log_file=mafft_log_file,
        verbose=verbose
    )
    
    if not align_result["success"]:
        return {
            "success": False,
            "error": f"MAFFT alignment failed: {align_result.get('error', 'Unknown error')}"
        }
    
    # Step 2: ClipKIT trimming
    if verbose:
        print(f"Trimming alignment with ClipKIT...")
    
    trim_result = run_clipkit_trimming(
        input_alignment=aligned_file,
        output_alignment=trimmed_file,
        mode=clipkit_mode,
        log_file=clipkit_log_file,
        verbose=verbose
    )
    
    if not trim_result["success"]:
        return {
            "success": False,
            "aligned_file": aligned_file,
            "error": f"ClipKIT trimming failed: {trim_result.get('error', 'Unknown error')}"
        }
    
    result = {
        "success": True,
        "aligned_file": aligned_file,
        "trimmed_file": trimmed_file
    }
    
    if mafft_log_file and "log_file" in align_result:
        result["mafft_log"] = align_result["log_file"]
    if clipkit_log_file and "log_file" in trim_result:
        result["clipkit_log"] = trim_result["log_file"]
    
    return result

