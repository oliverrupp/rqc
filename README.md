# RQC — RNA-seq Quality Control Pipeline

RQC is a Snakemake-based RNA-seq quality control workflow that generates an interactive HTML report covering sequencing quality, library complexity, RNA degradation, replicate consistency, and differential expression diagnostics.

## Requirements

Only **Conda** (Miniconda or Anaconda) must be installed and the user must have the rights to create an environment.

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

---

## Documentation
- [Usage](docs/usage.md)
- [Input Data Format](docs/input-data.md)
- [Workflow](docs/workflow.md)
- [Using a HPC Cluster](docs/hpc.md)
- [QC Metrics](docs/qc-metrics.md)
- [Metric Interpretation Guide](docs/interpretation.md)
- [Benchmarking and Performance](docs/benchmarking.md)
- [FAQ](docs/faq.md)
