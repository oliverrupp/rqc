# Output Files

## Output Directory Structure

The pipeline produces the following output structure for each subproject:

```
subproject_name/
├── report/
│   └── samples.report.html              # Interactive HTML QC report
├── results/
│   ├── salmon/                          # Transcript quantification (tximeta)
│   │   └── {sample}/quant.sf            # Salmon quantification files
│   ├── salmon_rrna/                     # rRNA quantification
│   ├── salmon_quantiles/                # Quantile-stratified quantification
│   ├── trimmed/                         # Trimmed sequencing reads (FastP)
│   │   ├── {sample}.json                # FastP JSON report
│   │   ├── {sample}_1.fq.gz             # Trimmed forward reads
│   │   └── {sample}_2.fq.gz             # Trimmed reverse reads
│   ├── falco/                           # Sequencing QC reports (Falco)
│   │   ├── trimmed/
│   │   └── untrimmed/
│   └── index/                           # Salmon index and metadata
│       ├── salmon/                      # Salmon index
│       └── BFC/                         # BiocFileCache
├── benchmark/                           # Execution time/memory benchmarks
├── logs/                                # Snakemake log files
└── tsv/
    └── samples/
        ├── complexity.tsv               # Complexity metrics
        ├── coverage_skewness.tsv        # Degradation metrics (skewness)
        ├── gene_body_coverage.tsv       # Gene body coverage (10 quantiles)
        ├── pca.tsv                      # PCA scores (PC1-PC5)
        ├── pca_outlier.tsv              # Outlier detection results
        ├── condition_df.tsv             # Condition QC scores
        ├── gene_detection.tsv           # Gene detection saturation curve
        ├── size_factor_qc.tsv           # DESeq2 size factors
        ├── g_count_matrix.tsv           # Raw gene-level counts
        ├── g_TPM.tsv                    # Gene-level TPM
        ├── g_TMM.tsv                    # Gene-level TMM (EdgeR)
        ├── g_geTMM.tsv                  # Gene-level geTMM
        ├── g_vst.tsv                    # Gene-level VST (DESeq2)
        ├── g_dispersion.tsv             # Gene dispersion estimates
        ├── t_count_matrix.tsv           # Raw transcript-level counts
        ├── t_TPM.tsv                    # Transcript-level TPM
        ├── t_TMM.tsv                    # Transcript-level TMM
        ├── t_geTMM.tsv                  # Transcript-level geTMM
        ├── t_vst.tsv                    # Transcript-level VST
        └── t_dispersion.tsv             # Transcript dispersion estimates
```

## Main Output: Interactive HTML Report

**Location**: `report/samples.report.html`

The interactive HTML report includes:

1. **Summary Tab**
   - Overview dashboard with traffic light metrics (Complexity, Degradation, Condition)
   - Per-sample QC scores

2. **Complexity Tab**
   - Read assignment stacked bar chart
   - Top N transcripts metrics table
   - Gene detection saturation curves

3. **Degradation Tab**
   - Gene body coverage heatmap (5' to 3' quantiles)
   - Coverage bias skewness metric

4. **PCA Tab**
   - Sample correlation heatmap
   - PCA scatter plots (PC1 vs PC2, etc.)
   - PCA pair plots (all PCs)
   - Outlier detection visualization

5. **DGE Tab**
   - Size factor distributions
   - Dispersion plot (gene-wise vs fitted)
   - Density plots by condition

6. **Sequencing QC Tab**
   - Per-base quality scores (before/after trimming)
   - Per-sequence quality distribution
   - GC content across positions
   - Insert size distribution (paired-end)
   - Adapter content detection

## File Formats

### TSV Format
- Tab-separated values
- First row: column names
- First column: row names (gene/sample IDs)
- Easy to import into R, Python, Excel

### JSON Format
- FastP reports in JSON format
- Machine-readable for automated analysis

### GZIP Compression
- All FASTQ files: `.fq.gz` or `.fastq.gz`
- Use `gunzip` or R/Python to decompress

## Accessing Output Programmatically

### R
```R
# Read TSV file
complexity <- read.table("tsv/samples/complexity.tsv", header=TRUE, row.names=1, sep="\t")

# Read expression matrix
counts <- read.table("tsv/samples/g_count_matrix.tsv", header=TRUE, row.names=1, sep="\t")
vst <- read.table("tsv/samples/g_vst.tsv", header=TRUE, row.names=1, sep="\t")
```

### Python
```python
import pandas as pd

# Read TSV file
complexity = pd.read_csv("tsv/samples/complexity.tsv", sep="\t", index_col=0)

# Read expression matrix
counts = pd.read_csv("tsv/samples/g_count_matrix.tsv", sep="\t", index_col=0)
vst = pd.read_csv("tsv/samples/g_vst.tsv", sep="\t", index_col=0)
```