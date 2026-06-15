# RQC — RNA-seq Quality Control Pipeline

RQC is a Snakemake-based RNA-seq quality control workflow that generates an interactive HTML report covering sequencing quality, library complexity, RNA degradation, replicate consistency, and differential expression diagnostics.

## Requirements

Only [**Conda**](https://www.anaconda.com/docs/getting-started/miniconda/main) (Miniconda or Anaconda) must be installed and the user must have the rights to create an environment.

All other software dependencies are installed automatically by the `rqc` script.

## Installation

Clone the repository:

```bash
git clone <repository-url>
cd rqc
```

Make the wrapper executable:

```bash
chmod +x rqc
```

First time running the script will check for Conda and will install the main Conda environment.

```bash
./rqc --help
```


## Conda Environment

Required Conda environments will be install by snakemake locally in the `projects` folder.

The first execution may take several minutes while dependencies are installed.

Subsequent executions in the same project reuse existing environments.

---

## Documentation
- [Usage](docs/usage.md)
- [Input Data Format](docs/input-data.md)
- [Workflow/Output](docs/workflow.md)
- [Using a HPC Cluster](docs/hpc.md)
- [QC Metrics](docs/qc-metrics.md)
- [Metric Interpretation Guide](docs/interpretation.md)
- [Benchmarking](docs/benchmarking.md)
- [FAQ](docs/faq.md)

---

## Using GenXBrowser

**TODO**

---

## Help

```text
usage: rqc [-h] [-v] [-l]
           [-o ORGANISM] [--alignment] [--assembly]
           [--no-conda]
		   [--max-cpus MAX_CPUS] [--max-jobs MAX_JOBS] [--max-memory MAX_MEMORY]
		   [--hpc EXECUTOR] [--hpc-config HPC_CONFIG]
		   [-n] [-r] [-k] [--config CONFIG]
           [project_dir]

RNA Quality Control Pipeline

positional arguments:
  project_dir              Project directory (default: current working directory)

options:
  -h, --help               show this help message and exit
  -v, --validate           Validate the project folder and all subfolders
  -l, --list-organisms     List all valid organisms and exit
                           
  -o, --organism ORGANISM  Comma-separated list of organisms to run (default: all)
                           
  --alignment              Compute quantification on genome alignments
                              (default: use pseudo-alignments)
  --assembly               Run de-novo genome-guided assembly
                           Extends reference annotation if available
                              (default: user provided GTF file)
                           
  --no-conda               Do not use conda environments (default: use conda)
                           
  --max-cpus MAX_CPUS      Maximum CPUs for local execution (default: 8)
  --max-jobs MAX_JOBS      Maximum parallel jobs for HPC execution (default: 100)
  --max-memory MAX_MEMORY  Maximum available memory (default: 16 Gb)
  --hpc EXECUTOR           Run pipeline on HPC cluster (slurm, lsf, pbs). (default: local)
  --hpc-config HPC_CONFIG  Snakemake HPC profile/config file (YAML format)
                           
  -n, --dry-run            Perform a dry run without executing jobs
  -r, --rerun-incomplete   Rerun incomplete jobs
  -k, --keep-going         Keep going if jobs fail
  --config CONFIG          Snakemake config file (YAML format)

```
