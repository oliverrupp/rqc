#!/usr/bin/env python3
"""
RNA Quality Control (RQC)

A Python wrapper script for running the Snakemake-based RQC pipeline.
Handles input validation, configuration management, and Snakemake execution.

Usage:
    rqc.py [OPTIONS]

Example:
    rqc.py --no-conda --max-cpus 16
    rqc.py --hpc slurm --max-jobs 10 --hpc-config hpc_config.yaml
    rqc.py --list-organisms
    rqc.py /path/to/project -no-conda
"""

import argparse
import sys
import psutil
import subprocess
import logging
from pathlib import Path
from typing import List, Optional, Set


SALMON_MEMORY_USAGE_FACTOR = 13


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)-8s - %(message)s',
    datefmt="%Y-%m-%d %H:%M:%S"
)
logger = logging.getLogger(__name__)

logging.addLevelName( logging.ERROR, "\033[1;31m%-8s\033[1;0m" % logging.getLevelName(logging.ERROR))
logging.addLevelName( logging.WARNING, "\033[1;31m%-8s\033[1;0m" % logging.getLevelName(logging.WARNING))

class RQCValidator:
    """Validates project directory structure and input files."""
    
    REQUIRED_GENOME_FILE = "reference/genome.fa"
    REQUIRED_ANNOTATION_GTF_FILE = "reference/annotation.gtf"
    REQUIRED_ANNOTATION_GFF3_FILE = "reference/annotation.gff3"
    REQUIRED_SAMPLES_FILE = "reference/samples.tsv"
    READS_PATTERNS = [
        "reads/*_1.fq.gz",
        "reads/*_2.fq.gz",
        "reads/*_s.fq.gz"
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

    
    def find_organisms(self, verbose: bool = True) -> Set[str]:
        """Find all valid organism directories."""
        organisms = set()
        
        for item in self.project_dir.iterdir():
            if item.is_dir() and not item.name.startswith('.'):
                organism_path = item
                if verbose:
                    print("")
                    logger.info(f"Checking \'{organism_path.name}\'")
                
                # Check if it has the required structure
                if self._validate_organism(organism_path, verbose=verbose):
                    organisms.add(item.name)
                    if verbose:
                        logger.info("valid organism")
        if verbose:
            print("")
        
        if not organisms:
            if verbose:
                logger.error(f"No valid organisms found in {self.project_dir}")
                logger.error(f"Organisms must contain: {self.REQUIRED_GENOME_FILE}, "
                             f"{self.REQUIRED_ANNOTATION_GFF3_FILE} or {self.REQUIRED_ANNOTATION_GFF3_FILE}, {self.REQUIRED_SAMPLES_FILE}, "
                             f"and reads files")
            return set()

        if verbose:
            logger.info(f"Found {len(organisms)} valid organism(s): {sorted(organisms)}")
            
        return organisms
    
    def _validate_organism(self, organism_path: Path, verbose: bool = True) -> bool:
        """Validate a single organism has required files."""
        # Check required reference files
        genome_file = organism_path / self.REQUIRED_GENOME_FILE
        annotation_gtf_file = organism_path / self.REQUIRED_ANNOTATION_GTF_FILE
        annotation_gff3_file = organism_path / self.REQUIRED_ANNOTATION_GFF3_FILE
        reads_dir = organism_path / "reads"
        
        if not genome_file.exists():
            if verbose:
                logger.error(f"Missing {self.REQUIRED_GENOME_FILE}")
            return False
        
        if not (annotation_gff3_file.exists() or annotation_gtf_file.exists()):
            if verbose:
                logger.error(f"Missing {self.REQUIRED_ANNOTATION_GFF3_FILE} and {self.REQUIRED_ANNOTATION_GTF_FILE}")
            return False
        
        if not reads_dir.exists() or not reads_dir.is_dir():
            if verbose:
                logger.error("Missing reads directory")
            return False

        sample_files = self.get_sample_files(organism_path, verbose)
        
        if not sample_files:
            if verbose:
                logger.error("No valid samples*.tsv files found")
            return False
        
        # Check for at least one reads file
        reads_files = list(reads_dir.glob("*.fq.gz"))
        if not reads_files:
            if verbose:
                logger.error("No .fq.gz files found")
            return False
                
        return True
    
    def get_sample_names(self, organism: str) -> List[str]:
        """Return sample TSV filenames for an organism."""
        return sorted(f.name for f in self.get_sample_files(self.project_dir / organism)
    )

    def get_sample_files(self, organism_path: Path, verbose: bool = False) -> set[Path]:
        """Return all valid samples*.tsv files."""
        
        sample_files = sorted((organism_path / "reference").glob("samples*.tsv"))
        valid_samples = set()
        
        for sample_file in sample_files:
            if not self._validate_samples_tsv(sample_file, verbose):
                if verbose:
                    logger.warning(f"Excluding invalid sample file: {sample_file.name}")
            else:
                valid_samples.add(sample_file)

        if not valid_samples:
            logger.error("No invalid sample file found")

        return valid_samples

    def get_report_targets(self, organisms: Optional[List[str]] = None) -> List[str]:
        """Return report targets for all samples*.tsv files."""

        targets = []

        selected = (
            organisms
            if organisms is not None
            else sorted(self.find_organisms(verbose=False))
        )

        for organism in selected:
            ref_dir = self.project_dir / organism 

            for sample_file in self.get_sample_files(ref_dir):
                suffix = sample_file.stem[len("samples"):]

                if suffix:
                    targets.append(f"{organism}/report/samples{suffix}.report.html")
                else:
                    targets.append(f"{organism}/report/samples.report.html")

        return targets

    def get_organism_memory_requirements(self) -> dict:
        """Return memory requirements (GB) for all valid organisms."""
        result = {}

        for organism in sorted(self.find_organisms(verbose=False)):
            genome_file = (self.project_dir / organism / self.REQUIRED_GENOME_FILE)

            size_gb = genome_file.stat().st_size / (1024**3)
            required_mem_gb = size_gb * SALMON_MEMORY_USAGE_FACTOR

            result[organism] = {
                "genome_size_gb": size_gb,
                "required_mem_gb": required_mem_gb
            }

        return result

    @staticmethod
    def _validate_samples_tsv(samples_file: Path, verbose: bool = True) -> bool:
        """Validate samples.tsv has at least 'condition' and 'sample' columns."""
        try:
            with open(samples_file, 'r') as f:
                header = f.readline().strip().split('\t')
                if 'condition' not in header or 'sample' not in header:
                    if verbose:
                        logger.error(f"{samples_file.name} missing 'condition' or 'sample' column")
                    return False
            return True
        except Exception as e:
            logger.error(f"Error reading samples.tsv: {e}")
            return False


class RQCPipeline:
    """Manages RQC pipeline execution."""
    
    # Supported HPC executors in Snakemake 9
    SUPPORTED_EXECUTORS = {
        "slurm": "SLURM job scheduler",
        "lsf": "LSF job scheduler",
        "pbs": "PBS/Torque job scheduler",
        "slurm_singularity": "SLURM with Singularity",
        "lsf_singularity": "LSF with Singularity"
    }
    
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
    
    def validate_hpc_executor(self, executor: str) -> bool:
        """Validate that the specified HPC executor is supported."""
        if executor not in self.SUPPORTED_EXECUTORS:
            logger.error(f"Unsupported HPC executor: {executor}")
            logger.error(f"Supported executors: {', '.join(self.SUPPORTED_EXECUTORS.keys())}")
            return False
        
        logger.info(f"Using HPC executor: {executor} ({self.SUPPORTED_EXECUTORS[executor]})")
        return True
    
    @staticmethod
    def validate_hpc_config(hpc_config: Optional[Path]) -> bool:
        """Validate that HPC config file exists if provided."""
        if hpc_config:
            if not hpc_config.exists():
                logger.error(f"HPC config file not found: {hpc_config}")
                return False
            logger.info(f"Using HPC config file: {hpc_config}")
        
        return True
    
    def build_snakemake_command(
        self,
        execution_mode: str,
        use_conda: bool,
        max_cpus: int,
        max_memory: int,
        max_jobs: int,
        dry_run: bool,
        rerun_incomplete: bool,
        keep_going: bool,
        executor: Optional[str] = None,
        hpc_config: Optional[Path] = None,
        organisms: Optional[List[str]] = None,
        config_file: Optional[Path] = None
    ) -> List[str]:
        """Build the Snakemake command with all options (Snakemake 9 compatible)."""
        cmd = ["snakemake", "--benchmark-extended", "-s", str(self.snakemake_file)]
        
        # Add workflow directory
        cmd.extend(["--directory", str(self.project_dir)])

        if execution_mode != "hpc":
            mem = psutil.virtual_memory()
            mem_mb = mem.total // 1024**2

            if max_memory > mem_mb:
                logger.warning(f"requested max-memory {max_memory} exceeds total RAM {mem_mb}!")
                logger.warning(f"setting max-memory to {mem_mb}!")
                max_memory = mem_mb

            
        cmd.extend(["--resources", f"mem_mb={max_memory}"])
        cmd.extend(["--config", f"max_mem_mb={max_memory}"])
        
        # Add execution mode
        if execution_mode == "hpc":
            # Snakemake 9 uses --executor option
            cmd.extend(["--executor", str(executor)])
            
            # Add HPC config file (Snakemake profile)
            if hpc_config:
                cmd.extend(["--profile", str(hpc_config)])
            
            # Add max-jobs for HPC execution
            cmd.extend(["--max-jobs", str(max_jobs)])
        else:  # local (default)
            cmd.extend(["--cores", str(max_cpus)])
        
        # Add conda support
        if use_conda:
            cmd.append("--use-conda")
        
        # Add dry run
        if dry_run:
            cmd.append("--dry-run")

        # Add rerun incomplete
        if rerun_incomplete:
            cmd.append("--rerun-incomplete")

        # Add keep going
        if keep_going:
            cmd.append("--keep-going")
            
        # Add config file (in addition to profile for HPC)
        if config_file:
            cmd.extend(["--configfile", str(config_file)])

        # Add targets
        validator = RQCValidator(self.project_dir)
        targets = validator.get_report_targets(organisms)
        cmd.extend(targets)
        
        return cmd
    
    def run_snakemake(self, cmd: List[str]) -> int:
        """Execute the Snakemake command."""
        logger.info(f"Executing: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                cwd=str(self.project_dir)
            )
            return result.returncode
        except Exception as e:
            logger.error(f"Error running Snakemake: {e}")
            return 1


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="RNA Quality Control (RQC) Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  rqc.py --no-conda --max-cpus 16
  rqc.py --hpc slurm --max-jobs 100 --hpc-config profile.yaml
  rqc.py --hpc lsf --max-jobs 50 --organism subproj1,subproj2
  rqc.py --list-organisms
  rqc.py --dry-run
  rqc.py --validate
  rqc.py /path/to/project 

Supported HPC executors:
  slurm              SLURM job scheduler
  lsf                LSF job scheduler
  pbs                PBS/Torque job scheduler
  slurm_singularity  SLURM with Singularity
  lsf_singularity    LSF with Singularity
        """
    )
    
    parser.add_argument(
        "project_dir",
        nargs="?",
        default=None,
        help="Project directory (default: current working directory)"
    )
    
    # List organisms option
    parser.add_argument(
        "--list-organisms",
        action="store_true",
        help="List all valid organisms and exit"
    )
    
    # HPC execution (optional - defaults to local if not specified)
    parser.add_argument(
        "--hpc",
        type=str,
        metavar="EXECUTOR",
        help="Run pipeline on HPC cluster with specified executor (slurm, lsf, pbs). Default: local execution"
    )
    
    # Conda
    parser.add_argument(
        "--no-conda",
        action="store_true",
        help="Do not use conda environments (default: use conda)"
    )
    
    # Resource options
    parser.add_argument(
        "--max-cpus",
        type=int,
        default=8,
        help="Maximum CPUs for local execution (default: 8)"
    )

    # Resource options
    parser.add_argument(
        "--max-memory",
        type=int,
        default=16000,
        help="Maximum available memory (default: 16 Gb)"
    )
    
    parser.add_argument(
        "--max-jobs",
        type=int,
        default=100,
        help="Maximum parallel jobs for HPC execution (default: 100)"
    )
    
    # HPC-specific options
    parser.add_argument(
        "--hpc-config",
        type=Path,
        help="Snakemake HPC profile/config file (YAML format)"
    )
    
    # Organisms
    parser.add_argument(
        "--organism",
        type=str,
        help="Comma-separated list of organisms to run (default: all)"
    )
    
    # Dry run
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Perform a dry run without executing jobs"
    )

    # rerun incomplete
    parser.add_argument(
        "--rerun-incomplete",
        action="store_true",
        help="Rerun incomplete jobs"
    )

    # keep going
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Keep going if jobs fail"
    )

    # just validate
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Validate the project folder and all subfolders"
    )

    # Config file
    parser.add_argument(
        "--config",
        type=Path,
        help="Snakemake config file (YAML format)"
    )
    
    return parser.parse_args()


def list_organisms(project_dir: Path) -> int:
    """List all valid organisms in the project directory."""
    validator = RQCValidator(project_dir)

    if not validator.validate_project_structure():
        logger.error("Project structure validation failed")
        return 1

    all_organisms = validator.find_organisms()

    if not all_organisms:
        logger.error("No valid organisms found")
        return 1

    mem_info = validator.get_organism_memory_requirements()

    print(f"\nValid organisms in {project_dir}:")
    print("-" * 70)
    print(
        f"{'Organism':20s} "
        f"{'Genome (GB)':>12s} "
        f"{'Mem. Req. (GB)':>15s} "
        f"Sample files"
    )
    print("-" * 70)

    for organism in sorted(all_organisms):
        info = mem_info[organism]
        sample_files = validator.get_sample_names(organism)

        print(
            f"{organism:20s} "
            f"{info['genome_size_gb']:12.2f} "
            f"{info['required_mem_gb']:15.1f} "
            f"{sample_files[0]}"
        )

        for sample_file in sample_files[1:]:
            print(
                f"{'':20s} "
                f"{'':12s} "
                f"{'':15s} "
                f"{sample_file}"
            )
    
    max_required_mem_gb = max(x["required_mem_gb"] for x in mem_info.values())
    print(f'Maximum required memory: {max_required_mem_gb:5.2f} GB')
    print("-" * 70)

    return 0


def main():
    """Main entry point."""
    args = parse_arguments()

    # Get script directory
    script_dir = Path(__file__).parent
    
    # Determine project directory (use current working directory if not specified)
    project_dir = Path(args.project_dir) if args.project_dir else Path.cwd()
    
    # Handle list-organisms option
    if args.list_organisms:
        return list_organisms(project_dir)
    
    # Determine execution mode (defaults to local if --hpc not specified)
    execution_mode = "hpc" if args.hpc else "local"
    use_conda = not args.no_conda 
    
    logger.info("RNA Quality Control (RQC) Pipeline")
    logger.info(f"Project directory: {project_dir}")
    logger.info(f"Execution mode: {execution_mode}")
    logger.info(f"Conda support: {use_conda}")
    
    # Validate project structure
    logger.info("Validating project structure...")
    validator = RQCValidator(project_dir)
    
    if not validator.validate_project_structure():
        logger.error("Project structure validation failed")
        sys.exit(1)
    
    # Find and filter organisms
    all_organisms = validator.find_organisms()
    if not all_organisms:
        logger.error("No valid organisms found")
        sys.exit(1)
    
    # Parse organism filter if provided
    selected_organisms = None
    if args.organism:
        selected_organisms = [sp.strip() for sp in args.organism.split(",")]
        
        # Validate requested organisms exist
        invalid = set(selected_organisms) - all_organisms
        if invalid:
            logger.error(f"Invalid organism(s): {', '.join(invalid)}")
            logger.error(f"Available organisms: {', '.join(sorted(all_organisms))}")
            sys.exit(1)
        
        logger.info(f"Running selected organisms: {', '.join(selected_organisms)}")
    else:
        logger.info(f"Running all organisms: {', '.join(sorted(all_organisms))}")
    
    # Validate config file if provided
    if args.config:
        if not args.config.exists():
            logger.error(f"Config file not found: {args.config}")
            sys.exit(1)
        logger.info(f"Using config file: {args.config}")
    
    # Validate Snakemake file
    logger.info("Validating Snakemake configuration...")
    pipeline = RQCPipeline(project_dir, script_dir)
    
    if not pipeline.validate_snakemake_file():
        logger.error("Snakemake validation failed")
        sys.exit(1)
    
    # HPC-specific validation
    if execution_mode == "hpc":
        if not pipeline.validate_hpc_executor(args.hpc):
            logger.error("HPC executor validation failed")
            sys.exit(1)
        
        if not pipeline.validate_hpc_config(args.hpc_config):
            logger.error("HPC config validation failed")
            sys.exit(1)
    
    # Build Snakemake command
    logger.info("Building Snakemake command...")
    cmd = pipeline.build_snakemake_command(
        execution_mode=execution_mode,
        use_conda=use_conda,
        max_cpus=args.max_cpus,
        max_memory=args.max_memory,
        max_jobs=args.max_jobs,
        dry_run=args.dry_run,
        rerun_incomplete=args.rerun_incomplete,
        keep_going=args.keep_going,
        executor=args.hpc,
        hpc_config=args.hpc_config,
        organisms=selected_organisms,
        config_file=args.config
    )
    
    # Run Snakemake
    if args.dry_run:
        logger.info("Running in DRY RUN mode - no jobs will be executed")

    if args.validate:
        logger.info("validation OK")
        print("\nWill execute:")
        print(f"{" ".join(cmd)}")
        sys.exit(0)
        
    logger.info("Starting pipeline execution...")
    returnee = pipeline.run_snakemake(cmd)
    
    if returnee == 0:
        logger.info("Pipeline execution completed successfully")
    else:
        logger.error(f"Pipeline execution failed with return code {returnee}")
    
    sys.exit(returnee)


if __name__ == "__main__":
    main()
