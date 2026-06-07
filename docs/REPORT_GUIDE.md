# RQC Report Guide

Navigate and interpret the interactive HTML report.

## Report Overview

The RQC report is an interactive HTML5 dashboard organized into 7 main sections:

1. **Summary** - Quick QC overview and status
2. **Library** - Read assignment and complexity metrics
3. **RNA Degradation** - Gene body coverage and RNA quality
4. **Replicate Consistency** - PCA and correlation analysis
5. **DGE** - Differential expression diagnostics
6. **Sequencing QC** - Read quality and adapter analysis
7. **Benchmarks** - Computational resource usage
8. **Help** - Reference documentation

## Summary Section

### Overview
The first view when opening the report. Provides a color-coded summary matrix of sample quality.

### Key Metrics Displayed

| Metric | Column | Color Meaning |
|--------|--------|---------------|
| **Complexity** | Shows % reads in top 100 genes | Green = good, Red = poor |
| **Degradation** | Shows coverage skewness | Green = uniform, Red = skewed |
| **Condition** | Shows PCA clustering score | Green = tight clustering, Red = scattered |
| **Batch Factors** | Any additional metadata | Color-coded by factor value |

### How to Read

- **Green cells**: Metric is good
- **Yellow/orange cells**: Metric is borderline
- **Red cells**: Metric suggests quality issue

Hover over cells for exact values.

### Summary Description

The right panel explains the three main QC dimensions:
- What each metric measures
- How to interpret values
- What issues might cause poor metrics
- Where to find more details in the report

## Library Section

### Read Assignment

**Visualization:** Stacked bar chart showing read fate

**What it shows:**
- Vertical axis: Sample names (grouped by condition)
- Horizontal axis: Number of reads
- Colors: Different read categories
  - Green: Assigned reads (good)
  - Orange: Unassigned reads
  - Red: Low-quality/removed reads

**How to interpret:**
- **Tall green bars**: Most reads successfully quantified (good)
- **Short green bars with large gaps**: Poor mapping (investigate)
- **Large red sections**: Reads removed during trimming (check for quality issues)

**Actions:**
- Hover for exact read counts
- Click legend to toggle categories on/off
- Compare across conditions for consistency

### Library Metrics Table

Detailed metrics for each sample:

| Column | Meaning |
|--------|---------|
| `condition` | Experimental condition |
| `type` | Sequencing type (paired-end/single-end) |
| `orientation` | Library strandedness |
| `fragment_length` | Inferred insert size (bp) |
| `mapping_rate` | % reads mapped to transcriptome |
| `detected_genes` | Number of genes with ≥1 read |
| `top_100_transcripts` | % reads in top 100 genes (complexity) |
| `Duprate` | PCR duplication rate (%) |
| `normalized_AUC` | Library saturation metric |
| `shannon_entropy` | Gene distribution evenness |
| `effective_genes` | Exp(Shannon entropy) |

**How to use:**
- Sort by column (click header) to find outliers
- Compare within-condition for consistency
- Use `mapping_rate` to check annotation compatibility

### Gene Detection

**Visualization:** Line plot with multiple colored lines (one per sample)

**X-axis:** Fraction of reads (0-100%)
**Y-axis:** Number of genes detected

**What it shows:**
- How many genes are detected as sequencing depth increases
- Saturation point (plateau) indicates adequate sampling

**How to interpret:**
- **Steep rise early, plateau late**: Good coverage, could sequence more
- **Early plateau (< 50%)**: Saturation reached, may not need more reads
- **Spreading lines**: Samples saturating at different rates

**Grouping:** Facets group by "tail gain" (% gene gain from 75-100% reads)
- Samples with high tail gain benefit from more sequencing

## RNA Degradation Section

### Gene Body Coverage Heatmap

**Visualization:** Heatmap showing gene coverage profile

**X-axis:** Gene body position (10 quantiles: 5' → 3')
**Y-axis:** Sample names
**Color:** Coverage percentage (blue = low, red = high)

**How to interpret:**
- **Horizontal red line:** Uniform coverage (good)
- **3' bias (red at right):** Some 3' end bias (common, but problematic if strong)
- **5' bias (red at left):** Unusual - investigate
- **Patchy pattern:** May indicate artifacts

**Ordering:** Samples ordered by degradation skewness (worst first)

### Fragment Length Estimation

**Visualization:** Scatter plot showing inferred fragment lengths

**What it shows:**
- Estimated insert size for each sample
- Useful for verifying library characteristics

**Normal range:** 200-500 bp (typically ~250 bp peak)

**Outliers:** Investigate samples with very short (< 150 bp) or long (> 600 bp) fragments

## Replicate Consistency Section

### Sample Correlation Heatmap

**Visualization:** Heatmap of Pearson correlations between samples

**Top genes used:** Top 5000 most variable genes

**Color scale:** Blue = low correlation, Red = high correlation

**How to interpret:**
- **Red blocks along diagonal:** Strong within-condition correlations (good)
- **Blue blocks:** Poor correlation - investigate for outliers
- **Samples away from diagonal:** Potential batch effects

### Top 1000 Genes Similarity

**Visualization:** Heatmap of Jaccard similarity for top 1000 expressed genes

**What it shows:**
- Whether samples agree on which genes are highly expressed
- Sensitive to major expression changes

**Color scale:** Same as correlation heatmap

**Good sign:** Red blocks for within-condition samples

### PCA Plot (PC1 vs PC2)

**Visualization:** 2D scatter plot

**X-axis:** PC2 (2nd major variance component)
**Y-axis:** PC1 (major variance component)
**Colors:** Different colors for different conditions
**Point size:** Outlier influence (larger = more influential)

**How to interpret:**
- **Tight clustering by condition:** Good replicate consistency
- **Spread within condition:** Biological or technical variation
- **Large points:** Outlier samples (investigate)
- **Separation along axis:** Condition effect in that direction

**Interactive features:**
- Hover for sample name
- Can zoom/pan the plot

### PC Pairs Plot

**Visualization:** 6-panel grid showing PC1, PC2, PC3 pairwise relationships

**Diagonals:** Density distributions
**Off-diagonals:** Scatter plots

**Use:** Examine top 3 PCs to understand major variance structure

### PCA by Batch Factors

**Available if:** Batch factors in sample metadata

**Visualization:** Box plots of PC scores, grouped by batch factor

**What it shows:**
- PC score distribution for each batch factor level
- Helps identify confounding batch effects

**Good sign:** Similar distributions within same condition across batches
**Problem:** Different distributions for same condition across batches

## DGE Section

### Size Factor QC

**Visualization:** Density plot of log-fold-changes to geometric mean

**What it shows:**
- Distribution of normalization factors (size factors)
- Indicates whether samples have similar library compositions

**Good sign:** Density centered near 0 with similar spread
**Problem:** Bimodal distribution or wide spread

### Size Factor by Condition

**Visualization:** Faceted density plots (one per condition)

**Use:** Check for unusual composition effects specific to one condition

### Dispersion Estimates

**Visualization:** Scatter plot with trend line

**X-axis:** Mean normalized counts (log scale)
**Y-axis:** Gene-wise dispersion (log scale)

**Components:**
- **Black points:** Per-gene dispersion estimates
- **Red line:** Mean-variance relationship
- **Blue points:** Final fitted dispersions

**Interpretation:**
- **Smooth trend downward:** Normal behavior
- **Scatter around trend:** Expected (gene-to-gene variation)
- **Unusual patterns:** May indicate data issues

## Sequencing QC Section

### FastP Quality Curves

**Visualization:** Line plot showing quality score by position

**X-axis:** Position in read (bp)
**Y-axis:** Quality score (Phred)

**What it shows:**
- Quality score before and after trimming
- Helps identify problematic sequencing
- Facets show: untrimmed forward, untrimmed reverse, trimmed forward, trimmed reverse

**Normal pattern:**
- High quality at 5' end
- Gradual decline toward 3' end
- Quality improves after trimming

**Problem patterns:**
- Early drop: PCR quality issues
- Bimodal: Possible primer contamination

### FastP GC Content

**Visualization:** Line plot showing GC% by position

**What it shows:**
- GC content should be relatively stable across positions
- Deviations may indicate contamination or biases

**Expected:** Lines relatively flat

### FastP Insert Size Distribution

**Visualization:** Histogram of inferred insert sizes

**X-axis:** Insert size (bp)
**Y-axis:** Read count

**Use:** Verify that insert size matches library preparation

### Per-Base Sequence Quality (Falco)

**Visualization:** Line plot (if falco data available)

**Similar to FastP quality curves** but from falco report

### Per-Sequence Quality

**Visualization:** Density plot of per-read quality scores

**What it shows:**
- Distribution of average quality scores per sequence
- Good libraries show high quality across most reads

### Adapter Content

**Visualization:** Line plot showing adapter % by position

**Shows:** Detected adapters at each read position

**Good sign:** < 0.1% adapter content

**Problem:** Adapter content increasing toward 3' end

## Benchmarks Section

**Available if:** Benchmark files present

**Visualization:** Scatter plot showing computational resource usage

**X-axis:** Snakemake rule (job type)
**Y-axis:** Resource usage (memory or runtime)

**Facets:** Split by metric (memory vs runtime)

**Use:** Understand computational requirements for each processing step

## Help Section

Contains reference documentation within the report itself.

## Navigation Tips

### Moving Between Sections

- Click section name in left sidebar to navigate
- Report remembers scroll position

### Interactive Features

- **Hover:** Most plots show additional information on hover
- **Click legend:** Toggle data series on/off
- **Zoom/Pan:** Many plots allow zooming and panning
- **Download:** Most plots can be saved as PNG via plotly menu

### Exporting Data

- **TSV files:** All underlying data available as tab-separated files in `tsv/` directory
- **Plots:** Save as PNG/SVG using plot menu (camera icon)

### Reproducibility

- Report parameters stored in HTML (reproducible from same parameters)
- All source R code in `scripts/rqc.R` and `scripts/rqc.Rmd`
- Sample metadata embedded in report

## Common Workflows

### "Is my sample quality OK?"

1. Go to **Summary** tab
2. Look for red/orange in your sample's row
3. For each issue, click on relevant section:
   - Red **Complexity** → Go to **Library** tab
   - Red **Degradation** → Go to **RNA Degradation** tab
   - Red **Condition** → Go to **Replicate Consistency** tab

### "Are my replicates consistent?"

1. Go to **Replicate Consistency** tab
2. Look at correlation heatmap for within-condition blue blocks
3. Check PCA plot: do replicates cluster?
4. If issues, go to PCA by batch to check for confounding

### "Do I have contamination?"

1. Go to **Library** tab
2. Check **Read Assignment** chart for rRNA content
3. Check **Library Metrics** table for `mapping_rate` < 50%
4. If suspicious, check **Sequencing QC** for adapter content

### "Which samples are outliers?"

1. Go to **Replicate Consistency** tab
2. Look for large points in PCA plot (high outlier influence)
3. Check correlation heatmap for low correlations
4. Review in **Library** tab for anomalous metrics

## Troubleshooting Report Issues

**Report won't open:**
- Ensure modern browser (Chrome, Firefox, Safari, Edge)
- Check file is not corrupted (compare file size to expected)

**Plots not displaying:**
- Try refreshing browser
- Check browser console for errors
- Ensure JavaScript enabled

**Missing sections:**
- Some sections only appear if data available (e.g., Benchmarks if benchmark files exist)
- Check pipeline completed successfully

See [Troubleshooting Guide](TROUBLESHOOTING.md) for more help.
