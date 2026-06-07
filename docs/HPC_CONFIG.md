# HPC Configuration

## Overview

The RQC pipeline supports multiple HPC job schedulers through Snakemake 9's executor system. Configuration is managed via profile YAML files.

## Supported HPC Schedulers

| Executor | Scheduler | Typical Use |
|----------|-----------|-------------|
| `slurm` | SLURM | Most common in academic clusters |
| `lsf` | LSF/IBM | Enterprise environments |
| `pbs` | PBS/Torque | Some academic and HPC centers |
| `slurm_singularity` | SLURM + Singularity | Container-based workflows |
| `lsf_singularity` | LSF + Singularity | Container-based LSF systems |

## Basic HPC Execution

### SLURM

```bash
# Basic execution
rqc.py --hpc slurm --max-jobs 100

# With custom profile
rqc.py --hpc slurm --max-jobs 100 --hpc-config slurm_profile.yaml

# With specific subprojects
rqc.py --hpc slurm --max-jobs 50 --subproject exp1,exp2
```

### LSF

```bash
# Basic execution
rqc.py --hpc lsf --max-jobs 50

# With custom profile
rqc.py --hpc lsf --max-jobs 50 --hpc-config lsf_profile.yaml
```

### PBS/Torque

```bash
# Basic execution
rqc.py --hpc pbs --max-jobs 75

# With custom profile
rqc.py --hpc pbs --max-jobs 75 --hpc-config pbs_profile.yaml
```

## Profile Configuration

Create a YAML profile file to customize HPC execution.

### SLURM Profile Example

**slurm_profile.yaml**:
```yaml
executor: slurm
jobs: 100

default-resources:
  slurm_time: 120
  slurm_mem_mb: 4000
  slurm_partition: "normal"
  slurm_cpus: 4
  slurm_account: "myaccount"
```

### LSF Profile Example

**lsf_profile.yaml**:
```yaml
executor: lsf
jobs: 50

default-resources:
  mem_mb: 4000
  cpus: 4
  time: 120
  queue: "normal"
```

### PBS Profile Example

**pbs_profile.yaml**:
```yaml
executor: pbs
jobs: 75

default-resources:
  pbs_walltime: "02:00:00"
  pbs_mem_mb: 4000
  pbs_cpus: 4
  pbs_queue: "standard"
  pbs_account: "myproject"
```

## Monitoring Jobs

### SLURM
```bash
squeue -u $USER
squeue -u $USER -l
scancel JOBID
```

### LSF
```bash
bjobs
bjobs -l JOBID
bkill JOBID
```

### PBS
```bash
qstat -u $USER
qstat -f JOBID
qdel JOBID
```

## Best Practices

1. **Start small**: Test with a single subproject first
2. **Monitor resources**: Check `squeue` or `bjobs` while running
3. **Adjust incrementally**: Increase `--max-jobs` gradually
4. **Use dry-run**: Test configuration with `--dry-run`
5. **Log output**: Keep Snakemake logs for debugging