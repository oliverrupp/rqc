# Installation & Setup

## Requirements

Ensure you have the following dependencies installed:

- Python 3.7+
- Snakemake 9.0+
- Conda (for environment management)
- R 4.0+ with required packages

## R Package Dependencies

The pipeline requires the following R packages:
- tximeta
- edgeR
- DESeq2
- pracma
- jsonlite
- rmarkdown
- DT
- plotly
- heatmaply
- RColorBrewer
- GGally
- getopt
- gtools
- tidyverse
- filelock

## Installation

1. **Install Snakemake 9.0 or higher**

   ```bash
   conda install -c bioconda snakemake=9
   ```

2. **Install R packages** (optional, can be installed automatically via conda)

   ```R
   # In R console
   install.packages(c("tidyverse", "plotly", "heatmaply", "RColorBrewer", "GGally", "getopt", "gtools", "filelock"))
   
   # From BioConductor
   BiocManager::install(c("tximeta", "edgeR", "DESeq2", "pracma", "jsonlite", "rmarkdown", "DT"))
   ```

3. **Place rqc.py in your project**

   Place `rqc.py` in your project directory alongside the `rqc.smk` Snakemake workflow file.

4. **Verify installation**

   ```bash
   rqc.py --help
   ```

## Optional: Conda Environment

Create a dedicated conda environment for the RQC pipeline:

```bash
conda create -n rqc-pipeline python=3.10 snakemake=9
conda activate rqc-pipeline
```

Then install R packages into this environment:

```bash
conda install -c bioconda r-tidyverse r-plotly r-heatmaply r-rcolorbrewer r-ggally
conda install -c bioconda bioconductor-tximeta bioconductor-edger bioconductor-deseq2
```