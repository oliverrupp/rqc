# QC Metrics Explanation

In-depth guide to interpreting RQC quality metrics.

## Understanding the Three QC Dimensions

The RQC report summarizes quality using three main dimensions that capture different aspects of RNA-seq data quality:

### 1. Complexity

**What it measures:** How evenly reads are distributed across the transcriptome.

**Key metric:** Percentage of reads assigned to the top 100 genes

**Interpretation:**
- **Low complexity** (< 20%): Good - Reads spread across many genes
- **Moderate** (20-40%): Acceptable - Some genes more abundant
- **High complexity** (> 50%): Concerning - Few genes dominate

**What causes low complexity (issues)?**
1. **rRNA contamination** - Ribosomal RNA can comprise 50-90% of reads
2. **Poor RNA quality** - Degraded RNA biased toward stable transcripts
3. **Library construction issues** - PCR amplification bias or cDNA synthesis issues
4. **Tissue-specific expression** - Some tissues naturally have few highly expressed genes
5. **Single-strand bias** - Strand-specific issues in library preparation

**How to assess:**
1. Check rRNA contamination % in read mapping table
2. Look at gene detection curve - does saturation plateau?
3. Compare to expected gene count for organism
4. Review Shannon entropy and effective genes

### 2. Degradation

**What it measures:** RNA quality, specifically evidence of degradation.

**Key metric:** Gene body coverage skewness (5' to 3' bias)

**Interpretation:**
- **Near 0** (-0.5 to 0.5): Good - Uniform coverage across gene body
- **Positive (> 0.5)**: Concerning - 3' bias (reads overrepresented at 3' end)
- **Negative (< -0.5)**: Less common - 5' bias

**Why it matters:**
Degraded RNA is preferentially fragmented at random positions, which can create apparent 3' bias because:
- 3' end of remaining transcript can still be sequenced
- 5' end may be lost due to fragmentation
- Polyadenosine tail at 3' provides stability

**How to assess:**
1. Review gene body coverage heatmap - is coverage uniform?
2. Check fragment length estimates - unusually short?
3. Look at individual sample coverage profiles in RNA Degradation section
4. Examine Falco quality curves - any sudden quality drops?

**Expected patterns by organism:**
- **Eukaryotes**: Usually see some 3' bias due to poly-A capture
- **mRNA-Seq**: Can have stronger 3' bias than total RNA
- **dUTP/dTTP libs**: Often show more uniform coverage

### 3. Replicate Consistency

**What it measures:** How well biological replicates cluster together.

**Key metrics:**
- PCA sample clustering by condition
- Sample correlation heatmap
- Jaccard similarity of top 1000 genes
- Median distance to condition centroid

**Interpretation:**
- **Good (✓)**: Replicates cluster tightly by condition
- **Moderate (⚠)**: Some scatter within condition groups
- **Poor (✗)**: Samples don't cluster by condition, or batch effects dominate

**What causes clustering issues:**
1. **Biological outliers** - Sample contamination or mislabeling
2. **Technical issues** - Library prep failure, sequencing problems
3. **Batch effects** - Sample processing date, operator, batch effects
4. **Experimental design** - Unequal effect sizes between conditions

**How to assess:**
1. PCA plot: Do replicates of same condition cluster?
2. Correlation heatmap: Are within-condition correlations high (r > 0.8)?
3. PC by batch: Are batch factors orthogonal to condition?
4. Outlier flags: Are any samples flagged as outliers?

## Library Metrics in Detail

### Read Assignment Flow

Reads follow this path:
```
Raw Reads (100%)
    ↓
After FastP trimming (% low quality removed)
    ↓
Aligned to transcriptome via Salmon
    ├─ Unmapped (ambiguous or non-aligning)
    ├─ NoFeature (aligned but not in annotation)
    ├─ Assigned (successfully quantified)
    └─ rRNA (separately quantified against rRNA seqs)
```

**Interpretation:**
- **Assigned > 50%**: Generally good alignment
- **Assigned 30-50%**: May indicate annotation gaps or non-coding reads
- **Assigned < 30%**: Investigate: Check annotation coverage, look for rRNA contamination
- **rRNA > 10%**: Consider rRNA depletion for future experiments
- **Unmapped > 30%**: May indicate low quality reads or species mismatch

### Duplication Rate

**Definition:** Percentage of reads that appear to be exact PCR duplicates.

**Interpretation:**
- **< 5%**: Excellent - minimal PCR amplification bias
- **5-10%**: Good - acceptable duplication
- **10-30%**: Concerning - significant amplification bias
- **> 30%**: Problem - likely over-amplified library

**Causes:**
1. **Under-amplification** - Library not amplified enough (shouldn't happen)
2. **Over-amplification** - Too many PCR cycles
3. **Input level** - Very low RNA input increases duplication
4. **Library type** - Some prep methods intrinsically have more duplication

**Action:** If > 20%, consider reducing PCR cycles in future experiments.

### Fragment Length

**Definition:** Average inferred insert size (effective fragment length).

**Expected range:**
- **Typically 200-500 bp**
- **Peak around 250 bp** for standard Illumina adapters

**Variation by preparation:**
- **Shorter fragments (< 200 bp)**: Degraded RNA or size-selected small RNAs
- **Longer fragments (> 500 bp)**: High-quality RNA, longer initial fragmentation

**Use case:**
- Used by Salmon for improved quantification accuracy
- Can indicate RNA quality issues if unexpectedly short

## Gene Detection & Saturation

### Gene Detection Curve

**X-axis:** Fraction of sequencing reads (0-100%)
**Y-axis:** Number of genes detected (≥ 1 read)

**Interpretation:**
- **Steep initial rise**: Many genes detected at low coverage
- **Plateau**: Indicates saturation - diminishing returns from additional reads
- **Late plateau** (> 80%): Well-covered, may want more reads
- **Early plateau** (< 50%): Approaching saturation, may not need more reads

**Use:** Compare curves across samples to assess sampling depth adequacy.

### Saturation Metric

**Definition:** Area under the gene detection curve (normalized).

**Interpretation:**
- **High AUC (> 0.8)**: Good sampling coverage
- **Moderate (0.6-0.8)**: Acceptable coverage
- **Low (< 0.6)**: May want additional sequencing

## Degradation Metrics

### Gene Body Coverage Profile

**X-axis:** Position along gene (5' → 3', divided into 10 quantiles)
**Y-axis:** Percentage of reads at each position

**Ideal profile:** Relatively flat (uniform coverage)

**Problem profiles:**
- **3' bias**: High coverage at 3' end
  - Suggests RNA degradation
  - Exacerbated in poly-A selected libraries
  - Loss of 5' end during fragmentation

- **5' bias**: High coverage at 5' end
  - Rare in typical workflows
  - May indicate specific library prep issues

### Coverage Skewness Score

**Definition:** Quantifies asymmetry of gene body coverage

**Formula:** (Position_Max - Position_Min) / 7 × |Max - Min|

**Interpretation:**
- **-0.5 to 0.5**: Good - minimal bias
- **0.5 to 1.0**: Some 3' bias - acceptable
- **> 1.0**: Strong 3' bias - potential degradation
- **< -0.5**: 5' bias - unusual, investigate

## Size Factor & Dispersion (DESeq2 Diagnostics)

### Size Factors

**What:** Normalization factors that adjust for sequencing depth differences

**Distribution:** Should be similar across samples (typically 0.5-2.0)

**Interpretation:**
- **Tight distribution**: Good - samples well-balanced
- **Wide distribution**: Some samples have different depth/composition
- **Outliers > 2 or < 0.5**: Investigate - may indicate issues

**Plot:** Density plot of log-ratios to geometric mean
- **Centered near 0**: Balanced library sizes
- **Bimodal or multimodal**: May indicate sample composition differences

### Dispersion Estimates

**What:** Gene-wise variance estimation (DESeq2 parameter)

**Interpretation:**
- **Smooth trend**: Good - dispersion trend as expected
- **High scatter**: Normal - biological variance varies by gene
- **Unusual patterns**: Investigate - may indicate issues with count data

**Components shown:**
- **Black points**: Per-gene estimates
- **Blue points**: Final fitted dispersions
- **Red line**: Mean-variance relationship

## PCA Analysis

### Principal Components

**PC1 vs PC2:** Usually capture most variance

**What to look for:**
- **Condition separation**: Are conditions separated along one PC?
- **Replicate clustering**: Do replicates cluster tightly?
- **Outliers**: Samples distant from group centroid
- **Batch effects**: Do samples cluster by batch factor instead of condition?

### PCA Outlier Detection

**Method:** Leave-one-out distance to condition centroid

**Influence metric:** How much centroid moves when sample is removed

**Interpretation:**
- **Low influence**: Sample similar to other condition replicates
- **High influence**: Sample diverges from condition pattern
- **Flagged (is_outlier = TRUE)**: Influence > median + 3×MAD

**Action:** Investigate outlier samples for technical issues, mislabeling, or true biological variability.

## Batch Effect Assessment

**X-axis:** Batch factor values (e.g., batch, time point, operator)
**Y-axis:** PC rotation values

**Good sign:** Similar distributions across batch factor levels within condition

**Problem:** Different distributions for the same condition across batch levels
- Suggests batch effect confounding condition

**Remediation:** Include batch factor in DESeq2 design formula:
```r
design = ~batch + condition
```

## Quality Score Interpretation

### Per-Base Quality (FastQC/Falco)

**Y-axis:** Phred quality score (higher = better)
**X-axis:** Position in read

**Interpretation:**
- **Q > 30**: High confidence base calls (< 0.1% error rate)
- **Q 20-30**: Medium quality
- **Q < 20**: Low confidence

**Expected patterns:**
- **Illumina**: Quality stable until end, then drops
- **Good trimming**: Low-quality bases removed from 3' end

### Per-Sequence Quality

**X-axis:** Sequence quality score
**Y-axis:** Number of sequences

**Expected:** Most sequences at Q > 30

**Problem:** Bimodal distribution or peak at low quality
- Suggests sequencing issues or primer contamination

## Adapter Content

**Purpose:** Detects remaining sequencing adapters

**Action:** If > 0.1% adapter content at any position, may need more aggressive trimming

## Correlation Heatmap

**Color scale:** Red (high correlation) to Blue (low correlation)

**Interpretation:**
- **Within-condition**: Should be red (high correlation, r > 0.85)
- **Between-condition**: Depends on biology (often intermediate)
- **Any blue blocks**: Investigate - potential problem samples

**Use:** Quick visual identification of outlier samples.

## Best Practices for Interpretation

1. **Context matters**: Consider organism, tissue, preparation method
2. **Compare to controls**: Look at historical data from your lab
3. **Don't over-interpret small differences**: Focus on major outliers
4. **Investigate outliers**: Don't automatically exclude - understand why
5. **Combine metrics**: One outlier metric may be acceptable; multiple = problem

See [Report Guide](REPORT_GUIDE.md) for navigating the interactive report.
