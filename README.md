# RNA Quality Control (RQC)

A comprehensive Snakemake-based pipeline for RNA-seq quality control and statistical analysis. The RQC pipeline performs systematic quality assessment of RNA-seq data including read quality trimming, transcript quantification, rRNA contamination detection, and multi-dimensional QC analysis.

## Quick Start

### Installation & Setup

1. Clone the repository:
```bash
git clone <repository-url>
cd rqc
```

2. Create the required project directory structure (see [Project Structure](#project-structure) below)

3. Run the pipeline:
```bash
python rqc.py
```

### Basic Usage

```bash
# Run with all detected organisms
python rqc.py /path/to/project

# Run specific organisms only
python rqc.py --organism organism1,organism2

# Dry run to see what will be executed
python rqc.py --dry-run

# List all valid organisms in a project
python rqc.py --list-organisms

# Run on HPC with SLURM
python rqc.py --hpc slurm --max-jobs 100 --hpc-config slurm_profile.yaml

# Validate project structure
python rqc.py --validate
```

## Project Structure

The pipeline expects a hierarchical directory structure with one or more "organisms" (species or experiments). Each organism must contain:

```
project_root/
├── organism1/
│   ├── reference/
│   │   ├── genome.fa              # Reference genome (required)
│   │   ├── annotation.gtf         # Gene annotation in GTF format (recommended)
│   │   ├── annotation.gff3        # OR gene annotation in GFF3 format (will be converted)
│   │   ├── samples.tsv            # Sample metadata (required, see format below)
│   │   ├── samples_group1.tsv     # Optional: additional sample groupings
│   │   └── samples_groupN.tsv     # Multiple sample files supported
│   └── reads/
│       ├── sample1_1.fq.gz        # Paired-end read 1 (or _R1)
│       ├── sample1_2.fq.gz        # Paired-end read 2 (or _R2)
│       ├── sample2_1.fq.gz
│       ├── sample2_2.fq.gz
│       └── sample3_s.fq.gz        # Single-end reads
│
├── organism2/
│   ├── reference/
│   │   └── ...
│   └── reads/
│       └── ...
```

### Sample Metadata Format

The `samples.tsv` file is tab-delimited with required columns:

```
sample          condition       [batch_factor_1]    [batch_factor_2]
sample1         treatment       batch1              replicate1
sample2         treatment       batch1              replicate2
sample3         control         batch1              replicate1
sample4         control         batch1              replicate2
```

**Required columns:**
- `sample`: Unique sample identifier (must match read filenames)
- `condition`: Experimental condition for differential expression analysis

**Optional columns:**
- Any additional columns are treated as batch factors and included in QC analysis

## Features

The RQC pipeline provides comprehensive RNA-seq quality control including:

1. **Read Quality Assessment**
   - FastP-based read trimming and filtering
   - Per-base and per-sequence quality metrics (FastQC/Falco)
   - Adapter content analysis

2. **Transcript Quantification**
   - Salmon-based pseudo-alignment and quantification
   - Support for both rRNA-depleted and total RNA libraries
   - Quantile-based 5'/3' bias detection

3. **Contamination Detection**
   - rRNA contamination quantification
   - Barrnap-based rRNA identification

4. **Quality Metrics**
   - Library complexity analysis (Shannon entropy, gene detection)
   - RNA degradation assessment (gene body coverage)
   - Sample correlation and outlier detection
   - PCA-based replicate consistency
   - Read saturation curves

5. **Interactive Reports**
   - HTML5-based interactive dashboard
   - Multiple visualization formats (heatmaps, plots, tables)
   - Comprehensive help and documentation

For detailed documentation on features, metrics, and interpretation, see:
- [Features & Metrics](docs/FEATURES.md)
- [QC Metrics Explained](docs/QC_METRICS.md)
- [Report Guide](docs/REPORT_GUIDE.md)

## Command-Line Options

### Basic Options

```bash
-h, --help                  Show help message
--list-organisms            List all valid organisms in the project
--validate                  Validate project structure without running
--dry-run                   Show what would be executed without running
```

### Execution Modes

```bash
# Local execution (default)
python rqc.py --max-cpus 16 --max-memory 16000 --no-conda

# HPC execution
python rqc.py --hpc slurm --max-jobs 100 --hpc-config profile.yaml
```

### Resource Control

```bash
--max-cpus N                Maximum CPUs for local execution (default: 8)
--max-memory N              Maximum memory in MB (default: 16000)
--max-jobs N                Maximum parallel jobs on HPC (default: 100)
```

### HPC Options

```bash
--hpc EXECUTOR              HPC executor: slurm, lsf, pbs, slurm_singularity, lsf_singularity
--hpc-config FILE.yaml      Snakemake profile/config file for HPC
```

### Filtering Options

```bash
--organism ORG1,ORG2        Process only specified organisms (comma-separated)
--config FILE.yaml          Snakemake configuration file
--rerun-incomplete          Rerun incomplete jobs
--keep-going                Continue on job failures
```

## Execution Modes

### Local Execution

Run on your local machine:

```bash
python rqc.py --max-cpus 16 --max-memory 32000
```

### HPC Execution

Submit jobs to a cluster scheduler:

```bash
python rqc.py --hpc slurm --max-jobs 100 --hpc-config slurm_profile.yaml
```

Supported executors: `slurm`, `lsf`, `pbs`, `slurm_singularity`, `lsf_singularity`

For detailed setup instructions, see [HPC Setup Guide](docs/HPC_SETUP.md)

## Configuration

### Conda Environments

The pipeline uses conda environments for all tools. To use pre-installed tools instead:

```bash
python rqc.py --no-conda
```

The conda environments are defined in the `envs/` directory:
- `envs/salmon.yaml` - Salmon pseudo-alignment
- `envs/fastp.yaml` - Read trimming
- `envs/samtools.yaml` - BAM file processing
- `envs/STAR.yaml` - Genomic alignment (optional)
- `envs/R.yaml` - R analysis environment
- `envs/gffread.yaml` - GTF/GFF3 processing
- `envs/falco.yaml` - FastQC-like quality control
- `envs/barrnap.yaml` - rRNA detection
- `envs/scallop.yaml` - Transcript assembly (optional)

### Configuration Files

Customize analysis parameters via YAML config file:

```bash
python rqc.py --config custom_config.yaml
```

See [Configuration Guide](docs/CONFIGURATION.md) for available parameters.

## Output

The pipeline generates comprehensive outputs organized by organism and sample:

```
project_root/organism/
├── results/
│   ├── trimmed/              # Trimmed reads
│   ├── salmon/               # Transcript quantification
│   ├── salmon_rrna/          # rRNA quantification
│   ├── salmon_quantiles/     # 5'/3' bias quantification
│   ├── falco/                # FastQC-equivalent reports
│   │   ├── untrimmed/
│   │   └── trimmed/
│   ├── bam/                  # Aligned reads (optional)
│   └── index/                # Salmon/STAR indices
├── report/
│   └── samples.report.html   # Interactive QC dashboard
├── tsv/
│   └── samples/              # Detailed metrics (TSV format)
└── benchmark/                # Execution benchmarks
```

For detailed output descriptions, see [Output Guide](docs/OUTPUT.md)

## Pipeline Steps

The pipeline includes the following major steps:

1. **Validation** - Verify project structure and input files
2. **Reference Processing** - Index genome, extract transcripts, identify rRNA
3. **Quality Control** - Trim reads, assess base quality
4. **Quantification** - Align and quantify transcripts
5. **Analysis** - Calculate QC metrics
6. **Reporting** - Generate interactive HTML report

For a detailed workflow diagram and step descriptions, see [Workflow Guide](docs/WORKFLOW.md)

## Scripts

The pipeline includes several utility scripts:

- **rqc.py** - Main pipeline wrapper (CLI interface)
- **rqc.smk** - Snakemake workflow definition
- **scripts/rqc.R** - R-based analysis and metrics computation
- **scripts/rqc.Rmd** - R Markdown template for HTML report generation
- **scripts/render_md.R** - Render report from command line
- **scripts/split10.py** - Split transcripts into quantiles for 5'/3' bias detection
- **scripts/barrnap2gtf.py** - Convert barrnap rRNA predictions to GTF format

For detailed documentation on each script, see [Scripts Reference](docs/SCRIPTS.md)

## Troubleshooting

Common issues and solutions:

- **Project validation fails**: Check that all required files exist and are in the correct locations. Run `python rqc.py --list-organisms` to verify detected organisms.
- **Missing conda environments**: Ensure you have conda installed and available. The pipeline will create environments on first run.
- **Memory errors**: Increase `--max-memory` or reduce the number of threads. Large genomes may require up to 300GB RAM for indexing.
- **HPC submission failures**: Verify HPC profile configuration and cluster connectivity.

See [Troubleshooting Guide](docs/TROUBLESHOOTING.md) for more details.

## Architecture

The pipeline consists of:

1. **Python Wrapper** (`rqc.py`) - CLI interface, project validation, Snakemake command orchestration
2. **Snakemake Workflow** (`rqc.smk`) - Workflow definition with rules for each processing step
3. **R Analysis** (`scripts/rqc.R`) - Statistical analysis and metric computation
4. **R Markdown** (`scripts/rqc.Rmd`) - Report generation and visualization

For detailed architecture documentation, see [Architecture Guide](docs/ARCHITECTURE.md)

## Performance Considerations

- **Memory requirements**: Depends on genome size. Salmon indexing typically requires ~13x genome size in RAM
- **Runtime**: 200M paired-end reads typically process in 4-8 hours on 16 cores
- **Storage**: Plan for ~50-100GB per organism (results, indices, temporary files)

See [Performance Guide](docs/PERFORMANCE.md) for optimization tips.

## References

Key publications and tools:

- Salmon: Patro et al., Nature Methods 2017
- DESeq2: Love et al., Genome Biology 2014
- edgeR: Robinson et al., Bioinformatics 2010
- FastQC/Falco: Andrews et al.

## Support & Documentation

- [Features & Metrics](docs/FEATURES.md)
- [QC Metrics Explained](docs/QC_METRICS.md)
- [Report Guide](docs/REPORT_GUIDE.md)
- [Configuration Guide](docs/CONFIGURATION.md)
- [Output Guide](docs/OUTPUT.md)
- [Workflow Guide](docs/WORKFLOW.md)
- [Scripts Reference](docs/SCRIPTS.md)
- [HPC Setup Guide](docs/HPC_SETUP.md)
- [Troubleshooting Guide](docs/TROUBLESHOOTING.md)
- [Architecture Guide](docs/ARCHITECTURE.md)
- [Performance Guide](docs/PERFORMANCE.md)

## License

[Your license here]

## Citation

If you use RQC in your research, please cite:

[Citation information]
