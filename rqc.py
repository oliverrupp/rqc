#!/usr/bin/env python3
"""
RNA Quality Control (RQC) Pipeline Wrapper

A Python wrapper script for running the Snakemake-based RQC pipeline.
Handles input validation, configuration management, and Snakemake execution.

Usage:
    rqc.py [OPTIONS] [PROJECT_DIRECTORY]

Example:
    rqc.py /path/to/project --local --conda yes --max-cpus 16
    rqc.py /path/to/project --hpc --max-nodes 10 --subproject subproj1,subproj2
    rqc.py /path/to/project --list-subprojects
"""

import argparse
import os
import sys
import subprocess
import logging
from pathlib import Path
from typing import List, Optional, Set


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class RQCValidator:
    """Validates project directory structure and input files."""
    
    REQUIRED_GENOME_FILE = "reference/genome.fa"
    REQUIRED_ANNOTATION_FILE = "reference/annotation.gff3"
    REQUIRED_SAMPLES_FILE = "reference/samples.tsv"
    READS_PATTERNS = [
        "reads/*_1.fastq.gz",
        "reads/*_2.fastq.gz",
        "reads/*_s.fastq.gz"
    ]
    
    def __init__(self, project_dir: Path):
        self.project_dir = Path(project_dir).resolve()
    
    def validate_project_structure(self) -> bool:
        """Validate main project directory exists."""
        if not self.project_dir.exists():
            logger.error(f"Project directory does not exist: {self.project_dir}")
            return False
        
        if not self.project_dir.is_dir():
            logger.error(f"Project path is not a directory: {self.project_dir}")
            return False
        
        logger.info(f"Project directory found: {self.project_dir}")
        return True
    
    def find_subprojects(self) -> Set[str]:
        """Find all valid subproject directories."""
        subprojects = set()
        
        for item in self.project_dir.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                subproject_path = item
                
                # Check if it has the required structure
                if self._validate_subproject(subproject_path):
                    subprojects.add(item.name)
        
        if not subprojects:
            logger.error(f"No valid subprojects found in {self.project_dir}")
            logger.error(f"Subprojects must contain: {self.REQUIRED_GENOME_FILE}, "
                        f"{self.REQUIRED_ANNOTATION_FILE}, {self.REQUIRED_SAMPLES_FILE}, "
                        f"and reads files")
            return set()
        
        logger.info(f"Found {len(subprojects)} valid subproject(s): {sorted(subprojects)}")
        return subprojects
    
    def _validate_subproject(self, subproject_path: Path) -> bool:
        """Validate a single subproject has required files."""
        # Check required reference files
        genome_file = subproject_path / self.REQUIRED_GENOME_FILE
        annotation_file = subproject_path / self.REQUIRED_ANNOTATION_FILE
        samples_file = subproject_path / self.REQUIRED_SAMPLES_FILE
        reads_dir = subproject_path / "reads"
        
        if not genome_file.exists():
            logger.debug(f"Missing {self.REQUIRED_GENOME_FILE} in {subproject_path.name}")
            return False
        
        if not annotation_file.exists():
            logger.debug(f"Missing {self.REQUIRED_ANNOTATION_FILE} in {subproject_path.name}")
            return False
        
        if not samples_file.exists():
            logger.debug(f"Missing {self.REQUIRED_SAMPLES_FILE} in {subproject_path.name}")
            return False
        
        if not reads_dir.exists() or not reads_dir.is_dir():
            logger.debug(f"Missing reads directory in {subproject_path.name}")
            return False
        
        # Check for at least one reads file
        reads_files = list(reads_dir.glob("*.fastq.gz"))
        if not reads_files:
            logger.debug(f"No .fastq.gz files found in {subproject_path.name}/reads")
            return False
        
        # Validate samples.tsv has required columns
        if not self._validate_samples_tsv(samples_file):
            logger.debug(f"Invalid samples.tsv in {subproject_path.name}")
            return False
        
        return True
    
    def _validate_samples_tsv(self, samples_file: Path) -> bool:
        """Validate samples.tsv has at least 'condition' and 'sample' columns."""
        try:
            with open(samples_file, 'r') as f:
                header = f.readline().strip().split('\t')
                if 'condition' not in header or 'sample' not in header:
                    logger.debug(f"samples.tsv missing 'condition' or 'sample' column")
                    return False
            return True
        except Exception as e:
            logger.debug(f"Error reading samples.tsv: {e}")
            return False


class RQCPipeline:
    """Manages RQC pipeline execution."""
    
    def __init__(self, project_dir: Path, script_dir: Path):
        self.project_dir = Path(project_dir).resolve()
        self.script_dir = Path(script_dir).resolve()
        self.snakemake_file = self.script_dir / "rqc.smk"
    
    def validate_snakemake_file(self) -> bool:
        """Validate that rqc.smk exists in the script directory."""
        if not self.snakemake_file.exists():
            logger.error(f"Snakemake file not found: {self.snakemake_file}")
            return False
        
        logger.info(f"Snakemake file found: {self.snakemake_file}")
        return True
    
    def build_snakemake_command(
        self,
        execution_mode: str,
        use_conda: bool,
        max_cpus: int,
        max_nodes: int,
        dry_run: bool,
        subprojects: Optional[List[str]] = None,
        config_file: Optional[Path] = None
    ) -> List[str]:
        """Build the Snakemake command with all options."""
        cmd = ["snakemake", "-s", str(self.snakemake_file)]
        
        # Add workflow directory
        cmd.extend(["--directory", str(self.project_dir)])
        
        # Add execution mode
        if execution_mode == "hpc":
            # Example for HPC cluster (adjust based on your cluster type)
            cmd.append("--cluster")
            cmd.append("sbatch")  # or qsub, bsub, etc.
            cmd.extend(["--jobs", str(max_nodes)])
        else:  # local
            cmd.extend(["--cores", str(max_cpus)])
        
        # Add conda support
        if use_conda:
            cmd.append("--use-conda")
        
        # Add dry run
        if dry_run:
            cmd.append("--dryrun")
        
        # Add config file
        if config_file:
            cmd.extend(["--configfile", str(config_file)])
        
        # Add target rule or subprojects
        if subprojects:
            # Create output targets for specified subprojects
            targets = [f"'{sp}/report/samples.report.html'" for sp in subprojects]
            cmd.append(" ".join(targets))
        else:
            # Default: run report rule
            cmd.append("report")
        
        return cmd
    
    def run_snakemake(self, cmd: List[str]) -> int:
        """Execute the Snakemake command."""
        logger.info(f"Executing: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                " ".join(cmd),
                shell=True,
                cwd=str(self.project_dir)
            )
            return result.returncode
        except Exception as e:
            logger.error(f"Error running Snakemake: {e}")
            return 1


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="RNA Quality Control (RQC) Pipeline Wrapper",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rqc.py /path/to/project --local --conda yes --max-cpus 16
  rqc.py /path/to/project --hpc --max-nodes 10 --subproject subproj1,subproj2
  rqc.py /path/to/project --list-subprojects
  rqc.py /path/to/project --dry-run --local
        """
    )
    
    parser.add_argument(
        "project_dir",
        nargs="?",
        default=os.getcwd(),
        help="Project directory (default: current directory)"
    )
    
    # List subprojects option
    parser.add_argument(
        "--list-subprojects",
        action="store_true",
        help="List all valid subprojects and exit"
    )
    
    # Execution mode (not required if --list-subprojects is used)
    mode_group = parser.add_mutually_exclusive_group(required=False)
    mode_group.add_argument(
        "--local",
        action="store_true",
        help="Run pipeline on local machine"
    )
    mode_group.add_argument(
        "--hpc",
        action="store_true",
        help="Run pipeline on HPC cluster"
    )
    
    # Conda
    parser.add_argument(
        "--conda",
        choices=["yes", "no"],
        default="yes",
        help="Use conda environments (default: yes)"
    )
    
    # Resource options
    parser.add_argument(
        "--max-cpus",
        type=int,
        default=8,
        help="Maximum CPUs per job for local execution (default: 8)"
    )
    
    parser.add_argument(
        "--max-nodes",
        type=int,
        default=1,
        help="Maximum number of HPC nodes (default: 1)"
    )
    
    # Subprojects
    parser.add_argument(
        "--subproject",
        type=str,
        help="Comma-separated list of subprojects to run (default: all)"
    )
    
    # Dry run
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a dry run without executing jobs"
    )
    
    # Config file
    parser.add_argument(
        "--config",
        type=Path,
        help="Snakemake config file (YAML format)"
    )
    
    # Help is automatic
    
    return parser.parse_args()


def list_subprojects(project_dir: Path) -> int:
    """List all valid subprojects in the project directory."""
    validator = RQCValidator(project_dir)
    
    if not validator.validate_project_structure():
        logger.error("Project structure validation failed")
        return 1
    
    all_subprojects = validator.find_subprojects()
    
    if not all_subprojects:
        logger.error("No valid subprojects found")
        return 1
    
    # Print subprojects to stdout
    print(f"\nValid subprojects in {project_dir}:")
    print("-" * 60)
    for subproject in sorted(all_subprojects):
        print(f"  - {subproject}")
    print("-" * 60)
    print(f"Total: {len(all_subprojects)} subproject(s)")
    
    return 0


def main():
    """Main entry point."""
    args = parse_arguments()
    
    # Get script directory
    script_dir = Path(__file__).parent
    
    # Handle list-subprojects option
    if args.list_subprojects:
        return list_subprojects(Path(args.project_dir))
    
    # Execution mode is now required only if not listing subprojects
    if not args.local and not args.hpc:
        logger.error("Either --local or --hpc must be specified")
        sys.exit(1)
    
    # Determine execution mode
    execution_mode = "hpc" if args.hpc else "local"
    use_conda = args.conda == "yes"
    
    logger.info(f"RNA Quality Control (RQC) Pipeline Wrapper")
    logger.info(f"Execution mode: {execution_mode}")
    logger.info(f"Conda support: {use_conda}")
    
    # Validate project structure
    logger.info("Validating project structure...")
    validator = RQCValidator(args.project_dir)
    
    if not validator.validate_project_structure():
        logger.error("Project structure validation failed")
        sys.exit(1)
    
    # Find and filter subprojects
    all_subprojects = validator.find_subprojects()
    if not all_subprojects:
        logger.error("No valid subprojects found")
        sys.exit(1)
    
    # Parse subproject filter if provided
    selected_subprojects = None
    if args.subproject:
        selected_subprojects = [sp.strip() for sp in args.subproject.split(",")]
        
        # Validate requested subprojects exist
        invalid = set(selected_subprojects) - all_subprojects
        if invalid:
            logger.error(f"Invalid subproject(s): {', '.join(invalid)}")
            logger.error(f"Available subprojects: {', '.join(sorted(all_subprojects))}")
            sys.exit(1)
        
        logger.info(f"Running selected subprojects: {', '.join(selected_subprojects)}")
    else:
        logger.info(f"Running all subprojects: {', '.join(sorted(all_subprojects))}")
    
    # Validate config file if provided
    if args.config:
        if not args.config.exists():
            logger.error(f"Config file not found: {args.config}")
            sys.exit(1)
        logger.info(f"Using config file: {args.config}")
    
    # Validate Snakemake file
    logger.info("Validating Snakemake configuration...")
    pipeline = RQCPipeline(args.project_dir, script_dir)
    
    if not pipeline.validate_snakemake_file():
        logger.error("Snakemake validation failed")
        sys.exit(1)
    
    # Build Snakemake command
    logger.info("Building Snakemake command...")
    cmd = pipeline.build_snakemake_command(
        execution_mode=execution_mode,
        use_conda=use_conda,
        max_cpus=args.max_cpus,
        max_nodes=args.max_nodes,
        dry_run=args.dry_run,
        subprojects=selected_subprojects,
        config_file=args.config
    )
    
    # Run Snakemake
    if args.dry_run:
        logger.info("Running in DRY RUN mode - no jobs will be executed")
    
    logger.info("Starting pipeline execution...")
    returncode = pipeline.run_snakemake(cmd)
    
    if returncode == 0:
        logger.info("Pipeline execution completed successfully")
    else:
        logger.error(f"Pipeline execution failed with return code {returncode}")
    
    sys.exit(returncode)


if __name__ == "__main__":
    main()
