# Usage Guide

## Basic Syntax

```bash
rqc.py [OPTIONS] [PROJECT_DIRECTORY]
```

## Running in Current Directory

If no project directory is specified, the script uses the current working directory:

```bash
rqc.py [OPTIONS]
```

## Common Use Cases

### List Available Subprojects

```bash
# In current directory
rqc.py --list-subprojects

# In specific directory
rqc.py /path/to/project --list-subprojects
```

### Local Execution (Default)

```bash
# Run all subprojects with 16 CPUs
rqc.py --max-cpus 16

# Run with default CPU setting (8)
rqc.py

# Run specific subprojects
rqc.py --subproject subproj1,subproj2 --max-cpus 16

# Run with custom Snakemake config
rqc.py --config config.yaml --max-cpus 16

# Perform a dry run
rqc.py --dry-run --max-cpus 16
```

### HPC Execution

```bash
# SLURM execution
rqc.py --hpc slurm --max-jobs 100

# LSF execution
rqc.py --hpc lsf --max-jobs 50

# PBS execution
rqc.py --hpc pbs --max-jobs 50
```

### HPC with Custom Profile

```bash
# SLURM with custom profile
rqc.py --hpc slurm --max-jobs 100 --hpc-config profile.yaml

# LSF with custom profile
rqc.py --hpc lsf --max-jobs 50 --hpc-config lsf_profile.yaml

# PBS with custom profile
rqc.py --hpc pbs --max-jobs 50 --hpc-config pbs_profile.yaml
```

### Advanced Examples

```bash
# Run specific subprojects on SLURM
rqc.py --hpc slurm --max-jobs 50 --subproject subproj1,subproj2

# Run with Singularity containers
rqc.py --hpc slurm_singularity --max-jobs 100

# Disable conda environments
rqc.py --conda no --max-cpus 16

# Dry run on HPC
rqc.py --hpc slurm --dry-run

# Run in specific directory
rqc.py /data/rna_seq_project --hpc slurm --max-jobs 100
```

## Output

The pipeline generates:

1. **Interactive HTML Report**: `report/samples.report.html`
2. **Result Files**: Quantification, normalized counts, QC metrics
3. **TSV Files**: Detailed metrics in `tsv/samples/` directory

See OUTPUT.md for detailed information.

## Monitoring Progress

The script provides detailed logging:

```bash
# Run with verbose output (captured in snakemake logs)
rqc.py --max-cpus 16
```

Check Snakemake's output for:
- Current rule being executed
- Number of jobs completed
- Estimated time remaining
- Any warnings or errors

## Stopping the Pipeline

Press `Ctrl+C` to stop the pipeline gracefully. Snakemake will:
- Stop queuing new jobs
- Wait for running jobs to complete
- Clean up temporary files
- Exit cleanly

You can safely restart the pipeline later - completed steps will be skipped.