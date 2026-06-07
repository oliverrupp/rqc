# QC Metrics

The RQC pipeline computes a comprehensive set of quality control metrics and generates an interactive HTML report with visualizations. The following metrics are available:

## Complexity Metrics

These metrics assess the complexity and diversity of the RNA-seq library:

| Metric | Description | Interpretation |
|--------|-------------|------------------|
| **Top 10 Transcripts** | Percentage of reads attributed to top 10 genes | <5% = good complexity; >20% = low complexity |
| **Top 50 Transcripts** | Percentage of reads attributed to top 50 genes | <30% = good; >60% = concerning |
| **Top 100 Transcripts** | Percentage of reads attributed to top 100 genes | <50% = good; >80% = low complexity |
| **Duplication Rate** | Fraction of duplicate reads in sequencing library | <20% = good; >50% = high duplicates (library/PCR bias) |
| **Shannon Entropy** | Measure of transcript distribution diversity | Higher = more uniform distribution (better complexity) |
| **Saturation AUC** | Area under curve for saturation analysis | Closer to 1.0 = approaching saturation |
| **Gene Detection** | Cumulative fraction of reads across sorted genes | Shows how many additional genes detected at higher read counts |

## Read Mapping Metrics

Quality metrics from read trimming and alignment:

| Metric | Description | Notes |
|--------|-------------|-------|
| **Total Reads** | Number of sequencing reads before filtering | Raw count from sequencer |
| **Low Quality Reads** | Reads filtered due to low quality scores | Per FastP settings (default Q15) |
| **Too Short Reads** | Reads removed during trimming for being below minimum length | Per FastP settings (default 15bp) |
| **Too Long Reads** | Reads removed for exceeding maximum length | Per FastP settings (default 1000bp) |
| **N-reads** | Reads with too many ambiguous bases | Per FastP settings |
| **Unmapped Reads** | Reads that failed to map to reference transcriptome | Indicates sequence quality/contamination |
| **No Feature Reads** | Mapped reads not assigned to any feature | Intergenic or ambiguous alignments |
| **Assigned Reads** | Successfully mapped and feature-assigned reads | Used for quantification |
| **rRNA Reads** | Proportion of reads mapping to ribosomal RNA | Indicates ribosomal depletion efficiency |

## Degradation Analysis

Metrics to detect RNA degradation:

| Metric | Description | Red Flag |
|--------|-------------|----------|
| **Gene Body Coverage** | Distribution of reads across 5' to 3' regions of genes (10 quantiles) | Skewed towards 3' = degradation |
| **Coverage Skewness** | Asymmetry measure indicating preferential 3' or 5' bias | High skewness (>1) = significant degradation |

**Interpretation**:
- Healthy RNA: Relatively uniform coverage along genes
- Degraded RNA: Biased towards 3' end or very 5'-biased

## PCA & Sample Quality

Principal Component Analysis metrics:

| Metric | Description | Use Case |
|--------|-------------|----------|
| **PCA Scores** | PC1-PC5 scores for all samples (up to 90% variance explained) | Identify major sources of variation |
| **Sample Correlation** | Pearson correlation matrix between all samples | Check biological replicates cluster together |
| **Outlier Detection** | Leave-one-out influence scores to identify outlier samples | Flag problematic samples |
| **Condition QC Score** | Median distance from condition centroid in PCA space | Assess within-group consistency |
| **PCA Variance** | Percentage of variance explained by each principal component | Understand data structure |

**Outlier Criteria**: 
- Samples with influence > (median + 3 × MAD) are flagged as potential outliers
- Review these samples - they may have technical issues

## Normalization & Differential Expression Preparation

Metrics used in normalization and DE analysis:

| Metric | Description | Method |
|--------|-------------|--------|
| **Raw Count Matrix** | Gene and transcript-level count matrices | From Salmon quantification |
| **TPM (Transcripts Per Million)** | Abundance measure normalized for transcript length and library size | From Salmon + DESeq2 |
| **TMM (Trimmed Mean of M-values)** | EdgeR normalization factor (robust to highly-expressed genes) | EdgeR weighted trimmed mean |
| **geTMM** | Gene-level TMM normalization with length correction | EdgeR + gene annotation |
| **VST (Variance Stabilizing Transform)** | DESeq2 normalized values for variance stabilization | DESeq2 transformation |
| **Size Factors** | DESeq2-estimated size factors for between-sample normalization | DESeq2 median ratio |
| **Dispersion Estimates** | Gene-wise and fitted dispersions for negative binomial model | DESeq2 parameter estimation |

**Usage**:
- Use TPM for abundance-based comparisons
- Use VST counts for visualization and clustering
- Use size factors for DE analysis with DESeq2

## Sequencing QC (FastP)

Quality metrics from raw sequencing data:

| Metric | Description |
|--------|-------------|
| **Per-base Quality** | Mean quality score at each position in reads (before/after trimming) |
| **Per-sequence Quality** | Distribution of average quality scores across all reads |
| **GC Content** | Percentage of G+C bases at each position |
| **Insert Size Distribution** | Distribution of insert sizes (for paired-end sequencing) |
| **Adapter Content** | Detection and quantification of adapter sequences |

**Quality Score Interpretation**:
- Q20 = 1% error rate
- Q30 = 0.1% error rate
- Typical good quality: median Q30+ across reads

## Summary Metrics (Report Dashboard)

The HTML report displays an interactive summary dashboard with overall quality scores:

| Metric | Calculation | Range | Color Coding |
|--------|------------|-------|---------------|
| **Complexity** | Based on top 100 transcripts percentage | 0-100% | Green: <50%; Yellow: 50-70%; Red: >70% |
| **Degradation** | Based on coverage skewness (normalized) | 0-2 | Green: <0.3; Yellow: 0.3-0.6; Red: >0.6 |
| **Condition** | Median PCA distance within condition | 0-100 | Green: <1.5; Yellow: 1.5-3; Red: >3 |

**Color Scheme**:
- **Green**: High quality - no action needed
- **Yellow**: Acceptable quality - review if concerns
- **Red**: Low quality - investigate potential issues