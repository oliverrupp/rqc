# RQC — RNA-seq Quality Control Pipeline

RQC is a Snakemake-based RNA-seq quality control workflow that generates an interactive HTML report covering sequencing quality, library complexity, RNA degradation, replicate consistency, and differential expression diagnostics.

The main entry point is the `rqc` shell script, which automatically creates and manages the required Conda environment and executes the pipeline.

## Requirements

Only **Conda** (Miniconda or Anaconda) must be installed.

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

## Usage

Run the pipeline from a project directory:

```bash
rqc /path/to/project
```

Validate the project structure:

```bash
rqc --validate
```

List detected organisms:

```bash
rqc --list-organisms
```

Run a subset of organisms:

```bash
rqc --organism Arabidopsis,Oryza
```

Dry run:

```bash
rqc --dry-run
```

## Input Directory Structure

Projects may contain one or more organisms/species.

```text
project/
├── Arabidopsis/
│   ├── reads/
│   │   ├── sample1_1.fq.gz
│   │   ├── sample1_2.fq.gz
│   │   ├── sample2_1.fq.gz
│   │   └── sample2_2.fq.gz
│   │
│   └── reference/
│       ├── genome.fa
│       ├── annotation.gtf      # or GFF3, optional if assembly mode is used
│       ├── samples.tsv
│       └── samples_batch.tsv   # optional additional sample sheets
│
└── Oryza/
    └── ...
```

### Sample Table

At minimum, every `samples*.tsv` file must contain:

| column    | description                                |
| --------- | ------------------------------------------ |
| sample    | sample identifier matching FASTQ filenames |
| condition | biological condition                       |

Additional metadata columns are allowed and are used for batch-effect assessment.

Example:

```text
sample      condition   batch
sample1     control     A
sample2     control     B
sample3     treated     A
sample4     treated     B
```

## Output

For every sample table (`samples*.tsv`), the pipeline generates:

```text
organism/
└── report/
    └── samples.report.html
```

Interactive reports contain plots, tables, PCA visualizations, QC summaries, and sequencing diagnostics.

---

# QC Metrics

## Library Complexity

### Read Assignment

Distribution of reads across assigned features, rRNA, unmapped reads, low-quality reads, and filtering categories.

### Mapping Rate

Percentage of reads assigned to annotated features.

### Detected Genes

Number of genes/transcripts detected in a sample.

### Gene Detection Curve

Number of detected transcripts as a function of sequencing depth.

### Normalized AUC

Area under the gene-detection curve. Higher values indicate more efficient transcript discovery at lower sequencing depth.

### Tail Gain

Fraction of additional genes detected after ~75% of reads have been sampled. High values suggest sequencing depth has not yet saturated.

### Top-100 Transcript Fraction

Fraction of reads assigned to the 100 most abundant transcripts. Elevated values may indicate low library complexity.

### Duplication Rate

Estimated read duplication level.

### Shannon Entropy

Diversity of transcript abundance distribution.

### Effective Gene Count

Exponentiated Shannon entropy; interpretable as the effective number of expressed genes.

### Fragment Length

Estimated library fragment size.

---

## RNA Degradation

### Gene Body Coverage

Coverage distribution across transcript bodies.

### Degradation Score

5′/3′ coverage bias used to identify RNA degradation. Values near zero indicate uniform coverage.

---

## Replicate Consistency

### Sample Correlation

Pearson correlation between samples using highly variable genes.

### PCA

Principal component analysis of gene expression profiles to identify outliers and batch effects.

### PCA Distance to Centroid

Within-condition dispersion used to quantify replicate consistency.

### Top-1000 Gene Jaccard Similarity

Similarity of highly expressed genes between samples.

---

## Differential Expression Diagnostics

### Size Factors

Normalization factors used for count scaling.

### Dispersion Estimates

Gene-wise and fitted dispersion estimates used by differential expression methods.

---

## Sequencing QC

### Per-Base Quality

Base quality scores before and after trimming.

### Per-Sequence Quality

Distribution of read quality scores.

### GC Content

GC-content profiles across read positions.

### Adapter Content

Residual adapter contamination after trimming.

### Insert Size Distribution

Estimated fragment size distribution for paired-end libraries.

---

## Performance Metrics

Optional benchmark statistics collected from Snakemake:

* Runtime per workflow step
* Memory consumption per workflow step
* CPU utilization information
 
