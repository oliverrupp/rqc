# RNA Quality Control (RQC)

A Python command-line tool that orchestrates the execution of the Snakemake-based RNA Quality Control (RQC) pipeline. This wrapper handles project validation, configuration management, and simplifies pipeline execution on both local machines and HPC clusters.

## Quick Start

```bash
# List available subprojects
rqc.py --list-subprojects

# Run QC analysis locally
rqc.py --max-cpus 16

# Run on HPC cluster
rqc.py --hpc slurm --max-jobs 100
```

## Features

- **Easy Pipeline Execution**: Simple interface to run RNA-seq quality control analysis
- **Flexible Execution Modes**: Support for local execution and multiple HPC schedulers
- **Project Validation**: Automatic validation of project structure and required files
- **Subproject Management**: Run specific subprojects or all subprojects at once
- **HPC Support**: Compatible with SLURM, LSF, PBS, and Singularity-based executors
- **Snakemake 9 Compatible**: Uses modern Snakemake 9 executor system with profiles
- **Conda Integration**: Automatic conda environment management
- **Comprehensive QC Metrics**: Extensive quality control analysis with interactive visualizations

## Documentation

- [Installation & Setup](docs/INSTALLATION.md) - Requirements and installation instructions
- [Project Structure](docs/PROJECT_STRUCTURE.md) - How to organize your data
- [Usage Guide](docs/USAGE.md) - Command-line usage and examples
- [Parameters Reference](docs/PARAMETERS.md) - Detailed parameter documentation
- [QC Metrics](docs/QC_METRICS.md) - Available quality control metrics
- [Output Files](docs/OUTPUT.md) - Output directory structure and generated files
- [HPC Configuration](docs/HPC_CONFIG.md) - Setting up HPC execution
- [Troubleshooting](docs/TROUBLESHOOTING.md) - Common issues and solutions

## Requirements

- Python 3.7 or higher
- Snakemake 9.0 or higher
- Conda (if using `--conda yes`)
- R and required packages (tximeta, DESeq2, edgeR, etc.)
- For HPC: Appropriate scheduler installed (SLURM, LSF, or PBS)

## License

Refer to the repository LICENSE file for licensing information.

## Contact

For issues or questions, please open an issue on the GitHub repository.
