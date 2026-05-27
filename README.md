# RNA Quality Control (RQC) Pipeline Wrapper

A Python command-line tool that orchestrates the execution of the Snakemake-based RNA Quality Control (RQC) pipeline. This wrapper handles project validation, configuration management, and simplifies pipeline execution on both local machines and HPC clusters.

## Features

- **Easy Pipeline Execution**: Simple interface to run RNA-seq quality control analysis
- **Flexible Execution Modes**: Support for local execution and multiple HPC schedulers
- **Project Validation**: Automatic validation of project structure and required files
- **Subproject Management**: Run specific subprojects or all subprojects at once
- **HPC Support**: Compatible with SLURM, LSF, PBS, and Singularity-based executors
- **Snakemake 9 Compatible**: Uses modern Snakemake 9 executor system with profiles
- **Conda Integration**: Automatic conda environment management
- **Dry-run Mode**: Preview pipeline execution without running jobs

## Installation

Ensure you have the following dependencies installed:

- Python 3.7+
- Snakemake 9.0+
- Conda (for environment management)

Place `rqc.py` in your project directory alongside the `rqc.smk` Snakemake workflow file.

## Project Structure

The tool expects the following directory structure:

```
project_dir/
├── subproject1/
│   ├── reference/
│   │   ├── genome.fa
│   │   ├── annotation.gtf
│   │   └── samples.tsv
│   └── reads/
│       ├── sample1_1.fq.gz
│       ├── sample1_2.fq.gz
│       └── ...
├── subproject2/
│   └── ...
└── rqc.py
```

### Required Files per Subproject

- **reference/genome.fa** - Reference genome in FASTA format
- **reference/annotation.gtf** - Genomic annotations in GTF format
- **reference/samples.tsv** - Sample metadata with at least `sample` and `condition` columns
- **reads/** - Directory containing gzipped FASTQ files (`*.fq.gz`)

## Usage

### Basic Syntax

```bash
rqc.py [OPTIONS] [PROJECT_DIRECTORY]
```

### Running in Current Directory

If no project directory is specified, the script uses the current working directory:

```bash
rqc.py [OPTIONS]
```

## Parameters

### Positional Arguments

| Parameter | Description | Required | Default |
|-----------|-------------|----------|---------|
| `PROJECT_DIRECTORY` | Path to the project directory | No | Current working directory |

### Global Options

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `--list-subprojects` | flag | List all valid subprojects in the project and exit | - |
| `--conda` | choice | Use conda environments (`yes` or `no`) | `yes` |
| `--dry-run` | flag | Perform a dry run without executing jobs | - |
| `--config` | file path | Snakemake configuration file (YAML format) | - |

### Execution Mode Options

#### Local Execution (Default)

Used by default if `--hpc` is not specified.

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `--max-cpus` | integer | Maximum CPUs for local execution | `8` |

#### HPC Execution

| Parameter | Type | Description | Required |
|-----------|------|-------------|----------|
| `--hpc` | string | HPC executor to use | Yes (if using HPC) |
| `--max-jobs` | integer | Maximum parallel jobs for HPC execution | `100` |
| `--hpc-config` | file path | Snakemake HPC profile/config file (YAML format) | No |

**Supported HPC Executors:**

- `slurm` - SLURM job scheduler
- `lsf` - LSF job scheduler
- `pbs` - PBS/Torque job scheduler
- `slurm_singularity` - SLURM with Singularity containers
- `lsf_singularity` - LSF with Singularity containers

### Subproject Selection

| Parameter | Type | Description | Default |
|-----------|------|-------------|---------|
| `--subproject` | string | Comma-separated list of subproject names to run | All subprojects |

## Examples

### List Available Subprojects

```bash
# In current directory
rqc.py --list-subprojects

# In specific directory
rqc.py /path/to/project --list-subprojects
```

### Local Execution

```bash
# Run all subprojects with 16 CPUs
rqc.py --max-cpus 16

# Run specific subprojects
rqc.py --subproject subproj1,subproj2 --max-cpus 16

# Dry run
rqc.py --dry-run --max-cpus 16
```

### HPC Execution with SLURM

```bash
# Basic SLURM execution
rqc.py --hpc slurm --max-jobs 100

# SLURM with custom profile
rqc.py --hpc slurm --max-jobs 100 --hpc-config profile.yaml

# SLURM with specific subprojects
rqc.py --hpc slurm --max-jobs 50 --subproject subproj1,subproj2

# SLURM dry run
rqc.py --hpc slurm --dry-run
```

### HPC Execution with LSF

```bash
# Basic LSF execution
rqc.py --hpc lsf --max-jobs 50

# LSF with custom profile
rqc.py --hpc lsf --max-jobs 50 --hpc-config lsf_profile.yaml
```

### With Custom Snakemake Configuration

```bash
# Local execution with custom config
rqc.py --config config.yaml --max-cpus 16

# HPC execution with custom config
rqc.py --hpc slurm --config config.yaml --hpc-config profile.yaml
```

### Conda Disable

```bash
# Run without conda environments (use system Python)
rqc.py --conda no --max-cpus 16
```

## HPC Profile Configuration

When using HPC execution, you can provide a custom Snakemake profile via `--hpc-config`. This YAML file should contain cluster-specific settings.

Example `profile.yaml` for SLURM:

```yaml
executor: slurm
jobs: 100
default-resources:
  slurm_partition: "normal"
  slurm_time: 120
  slurm_mem_mb: 4000
```

See [Snakemake 9 documentation](https://snakemake.readthedocs.io/en/stable/executing/cloud.html) for detailed profile configuration options.

## Output

The pipeline produces quality control reports for each subproject:

```
subproject_name/
└── report/
    └── samples.report.html
```

## Exit Codes

- `0` - Pipeline executed successfully
- `1` - Pipeline failed (validation error or execution error)

## Troubleshooting

### No Valid Subprojects Found

Ensure each subproject directory contains all required files:
- `reference/genome.fa`
- `reference/annotation.gtf`
- `reference/samples.tsv`
- `reads/` directory with `.fq.gz` files
- `samples.tsv` with `sample` and `condition` columns

### Snakemake File Not Found

Ensure `rqc.smk` is located in the same directory as `rqc.py`.

### Invalid HPC Executor

Use one of the supported executors: `slurm`, `lsf`, `pbs`, `slurm_singularity`, or `lsf_singularity`.

### HPC Config File Not Found

Verify the path to your HPC profile file is correct and the file exists.

## Requirements

- Python 3.7 or higher
- Snakemake 9.0 or higher
- Conda (if using `--conda yes`)
- For HPC: Appropriate scheduler installed (SLURM, LSF, or PBS)
- For Singularity executors: Singularity/Apptainer installed

## License

Refer to the repository LICENSE file for licensing information.

## Contact

For issues or questions, please open an issue on the GitHub repository.
