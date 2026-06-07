# Parameters Reference

## Positional Arguments

| Parameter | Type | Description | Required | Default |
|-----------|------|-------------|----------|----------|
| `PROJECT_DIRECTORY` | string | Path to the project directory | No | Current working directory |

## Global Options

| Parameter | Type | Description | Default |
|-----------|------|-------------|----------|
| `--list-subprojects` | flag | List all valid subprojects in the project and exit | - |
| `--conda` | choice | Use conda environments (`yes` or `no`) | `yes` |
| `--dry-run` | flag | Perform a dry run without executing jobs | - |
| `--config` | file path | Snakemake configuration file (YAML format) | None |

### `--list-subprojects`
Lists all valid subprojects found in the project directory that meet the validation criteria.

**Usage**:
```bash
rqc.py --list-subprojects
rqc.py /path/to/project --list-subprojects
```

### `--conda`
Controls whether to use conda environments for dependencies.

**Values**:
- `yes` (default) - Use conda environments defined in the Snakefile
- `no` - Use system Python and installed packages

**Usage**:
```bash
rqc.py --conda yes --max-cpus 16
rqc.py --conda no --max-cpus 16
```

### `--dry-run`
Performs a dry run to preview what will be executed without actually running jobs.

**Usage**:
```bash
rqc.py --dry-run
rqc.py --hpc slurm --dry-run
```

### `--config`
Path to a custom Snakemake configuration file in YAML format.

**Usage**:
```bash
rqc.py --config config.yaml --max-cpus 16
rqc.py --hpc slurm --config config.yaml --max-jobs 100
```

## Execution Mode Options

### Local Execution (Default)

Used when `--hpc` is not specified.

| Parameter | Type | Description | Default |
|-----------|------|-------------|----------|
| `--max-cpus` | integer | Maximum CPUs for local execution | `8` |

#### `--max-cpus`
Maximum number of CPU cores to use for parallel execution.

**Usage**:
```bash
rqc.py --max-cpus 8
rqc.py --max-cpus 16
rqc.py --max-cpus 32
```

### HPC Execution

Used when `--hpc` is specified with an executor name.

| Parameter | Type | Description | Required | Default |
|-----------|------|-------------|----------|----------|
| `--hpc` | string | HPC executor to use | Yes (for HPC mode) | - |
| `--max-jobs` | integer | Maximum parallel jobs for HPC execution | No | `100` |
| `--hpc-config` | file path | Snakemake HPC profile/config file (YAML) | No | None |

#### `--hpc`
Specifies the HPC executor to use. This parameter triggers HPC execution mode.

**Supported Executors**:
- `slurm` - SLURM Workload Manager
- `lsf` - IBM Load Sharing Facility
- `pbs` - PBS/Torque Job Scheduler
- `slurm_singularity` - SLURM with Singularity container support
- `lsf_singularity` - LSF with Singularity container support

**Usage**:
```bash
rqc.py --hpc slurm --max-jobs 100
rqc.py --hpc lsf --max-jobs 50
rqc.py --hpc pbs --max-jobs 75
rqc.py --hpc slurm_singularity --max-jobs 100
```

#### `--max-jobs`
Maximum number of jobs that can run in parallel on the HPC cluster.

**Usage**:
```bash
rqc.py --hpc slurm --max-jobs 50
rqc.py --hpc lsf --max-jobs 100
```

#### `--hpc-config`
Path to a Snakemake profile configuration file for HPC-specific settings.

**Usage**:
```bash
rqc.py --hpc slurm --hpc-config profile.yaml
rqc.py --hpc lsf --hpc-config lsf_profile.yaml
```

## Subproject Selection

| Parameter | Type | Description | Default |
|-----------|------|-------------|----------|
| `--subproject` | string | Comma-separated list of subproject names to run | All subprojects |

### `--subproject`
Selects specific subprojects to run instead of running all available subprojects.

**Usage**:
```bash
# Run single subproject
rqc.py --subproject subproject1 --max-cpus 16

# Run multiple subprojects
rqc.py --subproject subproject1,subproject2 --max-cpus 16

# Run with HPC
rqc.py --hpc slurm --subproject subproject1,subproject2 --max-jobs 100

# List available subprojects first
rqc.py --list-subprojects
```

**Notes**:
- Subproject names are case-sensitive
- Names must match directory names in the project
- Use comma separation without spaces for multiple subprojects
- Invalid subproject names will cause an error with available options listed