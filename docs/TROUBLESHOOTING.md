# Troubleshooting

## Common Issues and Solutions

### Project Structure Issues

#### "No valid subprojects found"

**Cause**: Subproject directories don't meet validation requirements.

**Solution**:
1. Check each subproject contains all required files:
   ```bash
   ls -la subproject1/reference/
   # Should show: genome.fa, annotation.gtf, samples.tsv
   
   ls -la subproject1/reads/
   # Should show: *.fq.gz files
   ```

2. Verify samples.tsv has required columns:
   ```bash
   head subproject1/reference/samples.tsv
   # Must have 'sample' and 'condition' columns
   ```

3. Check file names match:
   - Sample names in `samples.tsv` should match read file names
   - Example: sample1 → sample1_1.fq.gz, sample1_2.fq.gz

#### "Snakemake file not found"

**Cause**: `rqc.smk` is not in the same directory as `rqc.py`.

**Solution**:
```bash
# Check files are in same directory
ls -la rqc.py rqc.smk

# If missing, copy from repository
cp /path/to/repo/rqc.smk ./
```

### Execution Issues

#### "Either --local or --hpc must be specified" (Old Error)

**Note**: This error no longer occurs - local is now the default.

**Current behavior**: 
- If no `--hpc` specified, runs locally
- Use `--hpc EXECUTOR` for HPC execution

#### Pipeline runs but produces no output

**Cause**: Check-only mode or interrupted execution.

**Solution**:
1. Verify pipeline actually ran:
   ```bash
   ls -la */report/
   # Should see *.report.html files
   ```

2. Check for errors in execution:
   ```bash
   # Look for error messages in console output
   rqc.py --max-cpus 16 2>&1 | tee run.log
   ```

3. Re-run if interrupted:
   ```bash
   rqc.py --max-cpus 16
   # Snakemake will resume from where it stopped
   ```

### HPC-Specific Issues

#### "Unsupported HPC executor"

**Cause**: Invalid executor name.

**Solution**:
```bash
# Use one of the supported executors:
rqc.py --hpc slurm --max-jobs 100      # ✓ Valid
rqc.py --hpc lsf --max-jobs 50         # ✓ Valid
rqc.py --hpc pbs --max-jobs 75         # ✓ Valid
rqc.py --hpc slurm_singularity ...     # ✓ Valid
# rqc.py --hpc pbs_singularity ...      # ✗ Not supported
```

#### "HPC config file not found"

**Cause**: Specified profile doesn't exist or path is wrong.

**Solution**:
```bash
# Create profile file
cat > slurm_profile.yaml << 'EOF'
executor: slurm
jobs: 100
default-resources:
  slurm_time: 120
  slurm_mem_mb: 4000
EOF

# Use with absolute or relative path
rqc.py --hpc slurm --hpc-config $(pwd)/slurm_profile.yaml
```

#### Jobs not submitting

**Cause**: Scheduler not available or credentials issue.

**Solution**:
```bash
# Test scheduler availability
squeue --version          # SLURM
bjobs                     # LSF
qstat                     # PBS

# Check cluster connectivity
ssh compute-node-1 echo "Connected"

# Verify account/project is correct in profile
cat slurm_profile.yaml | grep account
```

#### Jobs timing out

**Cause**: Time limit too short for job.

**Solution**:
```bash
# Increase time in profile
cat > slurm_profile.yaml << 'EOF'
executor: slurm
jobs: 100
default-resources:
  slurm_time: 240      # Increase from 120 minutes
  slurm_mem_mb: 4000
EOF

rqc.py --hpc slurm --hpc-config slurm_profile.yaml
```

### Conda Issues

#### "conda: command not found"

**Cause**: Conda not installed or not in PATH.

**Solution**:
```bash
# Verify conda is installed
which conda

# Initialize conda if needed
conda init

# Re-run rqc.py
rqc.py --max-cpus 16
```

#### "Failed to create conda environment"

**Cause**: Environment creation error (network, disk space, permissions).

**Solution**:
```bash
# Try without conda first
rqc.py --conda no --max-cpus 16

# Or clear conda cache
conda clean --all

# Or increase timeout
export CONDA_PKGS_DIRS=/tmp/conda_cache
rqc.py --max-cpus 16
```

### Memory Issues

#### "MemoryError" or "Out of memory"

**Cause**: Job exceeded allocated memory.

**Solution (Local)**:
```bash
# Reduce data or use fewer cores
rqc.py --max-cpus 4  # Down from 16
```

**Solution (HPC)**:
```yaml
# Increase memory in profile
default-resources:
  slurm_mem_mb: 8000  # Increase from 4000
```

#### Disk space errors

**Cause**: Insufficient disk space for outputs and temporary files.

**Solution**:
```bash
# Check disk usage
df -h /path/to/project

# Move to directory with more space
cd /large_storage/
rqc.py /path/to/project --max-cpus 16

# Or use alternate temp directory
export TMPDIR=/path/to/large/tmp
rqc.py --max-cpus 16
```

### R/Python Issues

#### "R packages not found"

**Cause**: R packages not installed in conda environment.

**Solution**:
```bash
# Install missing package
conda install -c bioconda r-tidyverse
conda install -c bioconda bioconductor-deseq2

# Or disable conda and use system R
rqc.py --conda no --max-cpus 16
```

#### Broken pipe or connection errors

**Cause**: Network issue or process crash during execution.

**Solution**:
```bash
# Re-run - Snakemake will resume from last completed step
rqc.py --max-cpus 16

# Or check for resource exhaustion
top
free -h
df -h
```

### File Permission Issues

#### "Permission denied" errors

**Cause**: Incorrect file permissions or directory ownership.

**Solution**:
```bash
# Check file permissions
ls -la reference/
ls -la reads/

# Fix permissions
chmod u+rw reference/samples.tsv
chmod u+rx reads/

# Or change directory ownership
chown -R $USER:$USER ./
```

### Diagnostic Steps

#### Enable verbose output

```bash
# Run with debug output
python -u rqc.py --max-cpus 16 2>&1 | tee debug.log

# Check log file for errors
tail -100 debug.log
```

#### Dry-run to validate

```bash
# Test without actual execution
rqc.py --dry-run --hpc slurm

# Shows what would run without submitting jobs
```

#### Validate project structure

```bash
# Check all files exist
for subproj in */; do
  echo "Checking $subproj"
  test -f "$subproj/reference/genome.fa" && echo "✓ genome.fa" || echo "✗ missing genome.fa"
  test -f "$subproj/reference/annotation.gtf" && echo "✓ annotation.gtf" || echo "✗ missing annotation.gtf"
  test -f "$subproj/reference/samples.tsv" && echo "✓ samples.tsv" || echo "✗ missing samples.tsv"
  test -d "$subproj/reads" && echo "✓ reads/" || echo "✗ missing reads/"
done
```

#### List available subprojects

```bash
rqc.py --list-subprojects

# This validates all subproject requirements
```

## Getting Help

1. **Check documentation**: Review relevant docs in `/docs` folder
2. **Run diagnostic**: Use steps above to identify the issue
3. **Check logs**: Review snakemake output and error messages
4. **Dry-run**: Validate configuration with `--dry-run`
5. **Ask for help**: Open an issue on GitHub with:
   - Error message
   - Command you ran
   - Project structure (output of `rqc.py --list-subprojects`)
   - System info (`python --version`, `snakemake --version`)
