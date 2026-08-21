import argparse
import shutil
import sys
import csv
import re
import subprocess
from typing import List
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments for ENA fastq downloading."""
    parser = argparse.ArgumentParser(
        description="Download fastq files from ENA using SRR IDs from a sample file",
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=35), 
        epilog="""
Examples:
  python download_ena.py samples.tsv
  python download_ena.py samples.tsv --tmpdir /custom/tmp
        """
    )
    
    parser.add_argument(
        "sample_file",
        nargs="?",
        default=None,
        help="TSV file containing sample data with 'sample' column containing SRR IDs"
    )
    
    parser.add_argument(
        "--tmpdir",
        default="/tmp",
        help="Temporary directory for fastq-dl to use (default: /tmp)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    return parser.parse_args()


def get_srr_ids(sample_file) -> list[str]:
    """
    Extract SRR IDs from a TSV sample file.
    
    Reads a TSV file and extracts all SRR IDs from the 'sample' column.
    Valid SRR IDs match the pattern [SED]RR[0-9]* (e.g., SRR123456, ERR789, DRR000).
    
    Args:
        sample_file: Path to TSV file with a 'sample' column
        
    Returns:
        List of valid SRR IDs found in the file
    """
    srr = re.compile('^[SED]RR[0-9]*$')

    ids = []
    
    with open(sample_file, newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            id = row['sample']
            m = srr.match(id)
            if m:
                ids.append(id)
        
    return ids


def download_srr(project_dir, id, tmpdir, verbose=False) -> int:
    """
    Download fastq files for a specific SRR ID using fastq-dl.
    
    Downloads files to a temporary directory and then moves them to the project's
    reads directory. Handles both single-end (SE) and paired-end (PE) reads.
    
    Files are renamed as:
    - SE reads: {id}_s.fq.gz
    - PE read 1: {id}_1.fq.gz
    - PE read 2: {id}_2.fq.gz
    - Metadata: {id}.metadata.tsv
    
    Args:
        project_dir: Root directory of the project
        id: SRR ID to download
        tmpdir: Temporary directory for fastq-dl to use
        verbose: If True, print status messages
        
    Returns:
        0 on success, 1 on failure
    """
    if verbose:
        print(f"Downloading SRR ID: {id}")
    
    cmd = ['fastq-dl', '-o', tmpdir, '-a', id]
    try:
        if verbose:
            print(f"  Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=not verbose, text=True)
        
        if result.returncode == 0:
            if verbose:
                print(f"  Download successful, processing files...")
            
            metadata = Path(tmpdir) / 'fastq-run-info.tsv'
            se_src = Path(tmpdir) / (id + '.fastq.gz')
            p1_src = Path(tmpdir) / (id + '_1.fastq.gz')
            p2_src = Path(tmpdir) / (id + '_2.fastq.gz')
            
            # Move metadata file if it exists
            if metadata.exists():
                dest = project_dir + '/reads/' + id + '.metadata.tsv'
                shutil.move(metadata, dest)
                if verbose:
                    print(f"  Moved metadata to {dest}")

            # Move paired-end reads if both exist
            if(p1_src.exists() and p2_src.exists()):
                dest1 = project_dir + '/reads/' + id + '_1.fq.gz'
                dest2 = project_dir + '/reads/' + id + '_2.fq.gz'
                shutil.move(p1_src, dest1)
                shutil.move(p2_src, dest2)
                if verbose:
                    print(f"  Moved paired-end reads to {dest1} and {dest2}")

            # Move single-end reads if they exist
            if(se_src.exists()):
                dest = project_dir + '/reads/' + id + '_s.fq.gz'
                shutil.move(se_src, dest)
                if verbose:
                    print(f"  Moved single-end reads to {dest}")
                
            if verbose:
                print(f"  Successfully processed {id}")
            return 0
        else:
            print(f"ERROR: fastq-dl failed for {id} with return code {result.returncode}")
            if verbose and result.stderr:
                print(f"  stderr: {result.stderr}")
            return 1
    except Exception as e:
        print("ERROR: " + str(e))
        return 1


def downloaded(project_dir, id) -> bool:
    """
    Check if fastq files for a given SRR ID have already been downloaded.
    
    Returns True if either single-end OR both paired-end reads exist in the project.
    
    Args:
        project_dir: Root directory of the project
        id: SRR ID to check
        
    Returns:
        True if downloaded, False otherwise
    """
    se_file = Path(project_dir + '/reads/' + id +'_s.fq.gz')
    p1_file = Path(project_dir + '/reads/' + id +'_1.fq.gz')
    p2_file = Path(project_dir + '/reads/' + id +'_2.fq.gz')

    return se_file.exists() or (p1_file.exists() and p2_file.exists())


def main():
    """Main entry point for the ENA download script."""
    args = parse_arguments()

    sample_file = args.sample_file
    tmpdir = args.tmpdir
    verbose = args.verbose

    if verbose:
        print(f"Starting ENA download script")
        print(f"  Sample file: {sample_file}")
        print(f"  Temp directory: {tmpdir}")

    # Extract project directory from sample file path
    # Assumes sample_file is at: project_dir/reference/something/samples.tsv
    project_dir = re.sub('/reference/.*', '', sample_file)
    
    if verbose:
        print(f"  Project directory: {project_dir}")

    # Get list of SRR IDs from sample file
    srr_ids = get_srr_ids(sample_file)
    
    if verbose:
        print(f"Found {len(srr_ids)} SRR IDs to process")

    # Download each SRR ID that hasn't been downloaded yet
    for id in srr_ids:
        if not downloaded(project_dir, id):
            download_srr(project_dir, id, tmpdir, verbose)
        elif verbose:
            print(f"Skipping {id} (already downloaded)")

    if verbose:
        print("Download script completed")

    
if __name__ == "__main__":
    main()
