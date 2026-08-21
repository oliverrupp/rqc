import argparse
import shutil
import sys
import csv
import re
import subprocess
from pathlib import Path


# Compile SRR ID regex pattern once at module level
SRR_PATTERN = re.compile(r'^[SED]RR\d+$')


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments for ENA fastq downloading."""
    parser = argparse.ArgumentParser(
        description="Download fastq files from ENA using SRR IDs from a sample file",
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=35), 
        epilog="""
Examples:
  python download_ena.py samples.tsv
  python download_ena.py samples.tsv --tmpdir /custom/tmp
  python download_ena.py samples.tsv --project-dir /path/to/project -v
        """
    )
    
    parser.add_argument(
        "sample_file",
        help="TSV file containing sample data with 'sample' column containing SRR IDs"
    )
    
    parser.add_argument(
        "--tmpdir",
        default="/tmp",
        help="Temporary directory for fastq-dl to use (default: /tmp)"
    )
    
    parser.add_argument(
        "--project-dir",
        default=None,
        help="Project directory (default: auto-detect from sample_file path)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose output"
    )

    return parser.parse_args()


def validate_inputs(sample_file: str, tmpdir: str, project_dir: str) -> None:
    """
    Validate that all required inputs are valid and accessible.
    
    Args:
        sample_file: Path to sample TSV file
        tmpdir: Temporary directory path
        project_dir: Project directory path
        
    Raises:
        FileNotFoundError: If sample_file or project_dir don't exist
        ValueError: If tmpdir is not writable or other validation fails
    """
    sample_path = Path(sample_file)
    if not sample_path.exists():
        raise FileNotFoundError(f"Sample file not found: {sample_file}")
    
    if not sample_path.is_file():
        raise ValueError(f"Sample file is not a regular file: {sample_file}")
    
    tmpdir_path = Path(tmpdir)
    if not tmpdir_path.exists():
        raise FileNotFoundError(f"Temp directory does not exist: {tmpdir}")
    
    if not os.access(tmpdir, os.W_OK):
        raise ValueError(f"Temp directory is not writable: {tmpdir}")
    
    project_path = Path(project_dir)
    if not project_path.exists():
        raise FileNotFoundError(f"Project directory does not exist: {project_dir}")
    
    reads_dir = project_path / 'reads'
    if not reads_dir.exists():
        try:
            reads_dir.mkdir(parents=True, exist_ok=True)
        except OSError as e:
            raise ValueError(f"Cannot create reads directory: {e}")


def extract_project_dir(sample_file: str) -> str:
    """
    Extract project directory from sample file path.
    
    Assumes sample_file is at: project_dir/reference/something/samples.tsv
    
    Args:
        sample_file: Path to sample file
        
    Returns:
        Extracted project directory
        
    Raises:
        ValueError: If path doesn't match expected structure
    """
    match = re.search(r'^(.+?)/reference/', sample_file)
    if not match:
        raise ValueError(
            f"Sample file path doesn't match expected structure "
            f"(project_dir/reference/.../samples.tsv): {sample_file}"
        )
    return match.group(1)


def get_srr_ids(sample_file: str) -> list[str]:
    """
    Extract SRR IDs from a TSV sample file.
    
    Reads a TSV file and extracts all SRR IDs from the 'sample' column.
    Valid SRR IDs match the pattern [SED]RR[0-9]+ (e.g., SRR123456, ERR789, DRR000).
    
    Args:
        sample_file: Path to TSV file with a 'sample' column
        
    Returns:
        List of valid SRR IDs found in the file
        
    Raises:
        FileNotFoundError: If sample_file doesn't exist
        ValueError: If 'sample' column is missing or file is invalid
    """
    try:
        with open(sample_file, newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            
            if reader.fieldnames is None or 'sample' not in reader.fieldnames:
                raise ValueError(
                    f"TSV file must contain a 'sample' column. "
                    f"Found columns: {reader.fieldnames}"
                )
            
            # Use list comprehension for cleaner filtering
            return [
                row['sample'] 
                for row in reader 
                if row['sample'] and SRR_PATTERN.match(row['sample'])
            ]
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Sample file not found: {sample_file}") from e
    except (csv.Error, KeyError) as e:
        raise ValueError(f"Error parsing sample file: {e}") from e


def download_srr(project_dir: str, srr_id: str, tmpdir: str, verbose: bool = False) -> bool:
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
        srr_id: SRR ID to download
        tmpdir: Temporary directory for fastq-dl to use
        verbose: If True, print status messages
        
    Returns:
        True on success, False on failure
    """
    if verbose:
        print(f"Downloading SRR ID: {srr_id}")
    
    cmd = ['fastq-dl', '-o', tmpdir, '-a', srr_id]
    try:
        if verbose:
            print(f"  Executing: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        
        if result.returncode != 0:
            print(f"ERROR: fastq-dl failed for {srr_id} with return code {result.returncode}")
            if result.stderr:
                print(f"  stderr: {result.stderr}")
            return False
        
        if verbose:
            print(f"  Download successful, processing files...")
        
        # Use Path for all file operations
        tmpdir_path = Path(tmpdir)
        reads_dir = Path(project_dir) / 'reads'
        
        metadata = tmpdir_path / 'fastq-run-info.tsv'
        se_src = tmpdir_path / f'{srr_id}.fastq.gz'
        p1_src = tmpdir_path / f'{srr_id}_1.fastq.gz'
        p2_src = tmpdir_path / f'{srr_id}_2.fastq.gz'
        
        # Move metadata file if it exists
        if metadata.exists():
            dest = reads_dir / f'{srr_id}.metadata.tsv'
            shutil.move(str(metadata), str(dest))
            if verbose:
                print(f"  Moved metadata to {dest}")

        # Move paired-end reads if both exist
        if p1_src.exists() and p2_src.exists():
            dest1 = reads_dir / f'{srr_id}_1.fq.gz'
            dest2 = reads_dir / f'{srr_id}_2.fq.gz'
            shutil.move(str(p1_src), str(dest1))
            shutil.move(str(p2_src), str(dest2))
            if verbose:
                print(f"  Moved paired-end reads to {dest1.name} and {dest2.name}")

        # Move single-end reads if they exist
        elif se_src.exists():
            dest = reads_dir / f'{srr_id}_s.fq.gz'
            shutil.move(str(se_src), str(dest))
            if verbose:
                print(f"  Moved single-end reads to {dest.name}")
        
        if verbose:
            print(f"  Successfully processed {srr_id}")
        return True
        
    except OSError as e:
        print(f"ERROR: File operation failed for {srr_id}: {e}")
        return False
    except subprocess.TimeoutExpired:
        print(f"ERROR: fastq-dl timed out for {srr_id}")
        return False


def is_downloaded(project_dir: str, srr_id: str) -> bool:
    """
    Check if fastq files for a given SRR ID have already been downloaded.
    
    Returns True if either single-end OR both paired-end reads exist in the project.
    
    Args:
        project_dir: Root directory of the project
        srr_id: SRR ID to check
        
    Returns:
        True if downloaded, False otherwise
    """
    reads_dir = Path(project_dir) / 'reads'
    se_file = reads_dir / f'{srr_id}_s.fq.gz'
    p1_file = reads_dir / f'{srr_id}_1.fq.gz'
    p2_file = reads_dir / f'{srr_id}_2.fq.gz'

    return se_file.exists() or (p1_file.exists() and p2_file.exists())


def main() -> int:
    """
    Main entry point for the ENA download script.
    
    Returns:
        0 on success, 1 if any downloads failed, 2 on validation error
    """
    try:
        args = parse_arguments()
        
        sample_file = args.sample_file
        tmpdir = args.tmpdir
        verbose = args.verbose
        
        # Determine project directory
        project_dir = args.project_dir
        if not project_dir:
            project_dir = extract_project_dir(sample_file)
        
        if verbose:
            print(f"Starting ENA download script")
            print(f"  Sample file: {sample_file}")
            print(f"  Temp directory: {tmpdir}")
            print(f"  Project directory: {project_dir}")
        
        # Validate inputs
        try:
            validate_inputs(sample_file, tmpdir, project_dir)
        except (FileNotFoundError, ValueError) as e:
            print(f"ERROR: {e}")
            return 2
        
        # Get list of SRR IDs from sample file
        try:
            srr_ids = get_srr_ids(sample_file)
        except (FileNotFoundError, ValueError) as e:
            print(f"ERROR: {e}")
            return 2
        
        if not srr_ids:
            print("WARNING: No valid SRR IDs found in sample file")
            return 0
        
        if verbose:
            print(f"Found {len(srr_ids)} SRR IDs to process")

        # Download each SRR ID that hasn't been downloaded yet
        failed_count = 0
        for srr_id in srr_ids:
            if is_downloaded(project_dir, srr_id):
                if verbose:
                    print(f"Skipping {srr_id} (already downloaded)")
            else:
                success = download_srr(project_dir, srr_id, tmpdir, verbose)
                if not success:
                    failed_count += 1

        if verbose:
            print("Download script completed")
        
        if failed_count > 0:
            print(f"WARNING: {failed_count}/{len(srr_ids)} downloads failed")
            return 1
        
        return 0

    except KeyboardInterrupt:
        print("\nInterrupted by user")
        return 130
    except Exception as e:
        print(f"FATAL ERROR: {e}")
        return 1

    
if __name__ == "__main__":
    sys.exit(main())
