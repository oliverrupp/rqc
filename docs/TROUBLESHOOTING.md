# Troubleshooting

## Common Issues and Solutions

### Project Structure Issues

#### "No valid subprojects found"

**Cause**: Subproject directories don't meet validation requirements.

**Solution**:
1. Check each subproject contains all required files:
   ```bash
   ls -la subproject1/reference/
   ls -la subproject1/reads/
   ```

2. Verify samples.tsv has required columns:
   ```bash
   head subproject1/reference/samples.tsv
   ```

#### "Snakemake file not found"

**Cause**: `rqc.smk` is not in the same directory as `rqc.py`.

**Solution**:
```bash
ls -la rqc.py rqc.smk
```

### Execution Issues

#### Pipeline runs but produces no output

**Cause**: Check-only mode or interrupted execution.

**Solution**:
```bash
ls -la */report/
rqc.py --max-cpus 16 2>&1 | tee run.log
```

### HPC-Specific Issues

#### "Unsupported HPC executor"

**Cause**: Invalid executor name.

**Solution**:
```bash
# Use one of the supported executors:
rqc.py --hpc slurm --max-jobs 100      # Valid
rqc.py --hpc lsf --max-jobs 50         # Valid
rqc.py --hpc pbs --max-jobs 75         # Valid
```

#### Jobs timing out

**Cause**: Time limit too short for job.

**Solution**:
```yaml
# Increase time in profile
default-resources:
  slurm_time: 240      # Increase from 120 minutes
```

### Memory Issues

#### "MemoryError" or "Out of memory"

**Cause**: Job exceeded allocated memory.

**Solution (Local)**:
```bash
rqc.py --max-cpus 4
```

**Solution (HPC)**:
```yaml
default-resources:
  slurm_mem_mb: 8000  # Increase from 4000
```

### Diagnostic Steps

#### Enable verbose output

```bash
python -u rqc.py --max-cpus 16 2>&1 | tee debug.log
```

#### Dry-run to validate

```bash
rqc.py --dry-run --hpc slurm
```

#### List available subprojects

```bash
rqc.py --list-subprojects
```