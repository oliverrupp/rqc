# QC Metrics
## Summary Page

#### **Summary Metrics**

Three major factors can affect the quality of RNA-seq data and the reliability of downstream analyses. This page summarizes each factor using a single metric. Additional details and diagnostic plots are provided on the corresponding pages for each factor.


##### **1. Complexity**

Library complexity describes how evenly sequencing reads are distributed across genes. If a large fraction of reads originates from only a small number of highly expressed genes, fewer reads remain available for the quantification of other genes. Consequently, fewer genes can be measured with sufficient precision for reliable differential expression analysis.

Library complexity is measured here as the percentage of reads assigned to the 100 most abundant genes. Higher values indicate lower complexity.

Low complexity can result from several factors, including overrepresentation of a small number of transcripts, rRNA contamination, poor RNA quality, or technical issues during library preparation.

The **Library** page provides additional metrics that help identify the underlying cause of low complexity.


##### **2. Degradation**

RNA degradation prior to library preparation can introduce systematic biases into transcript quantification. Degraded RNA molecules are fragmented before sequencing, reducing the likelihood that transcript sequences are represented uniformly. As a result, transcript abundance estimates may become biased and less comparable across samples.

In high-quality libraries, read coverage is expected to be approximately uniform across the gene body. To assess degradation, the pipeline calculates the gene body coverage bias. Negative values indicate a 5′ bias, positive values indicate a 3′ bias, and values close to zero are expected for libraries with little or no degradation.

The **RNA Degradation** page provides detailed gene body coverage profiles for each sample.


##### **3. Replicate Consistency**

Biological replicates are used to estimate variability within the same experimental condition and increase the reliability of differential expression analyses. Replicates belonging to the same condition are generally expected to exhibit similar gene expression patterns.

To assess replicate consistency, principal component analysis (PCA) is performed on the expression matrix. For each condition, the mean distance of samples to the centroid of their condition group is calculated. Larger values indicate greater variability among replicates and may suggest the presence of outlier samples, sample mix-ups, or other technical issues.

If additional metadata variables (batch factors) are available, the same metric is calculated for each factor. When samples cluster more strongly by a batch factor than by the biological condition of interest, this may indicate a batch effect that should be considered during downstream analyses.

The **PCA** page contains both the sample correlation heatmap and PCA plots, which can be used to identify outlier samples and potential batch effects.


---

## Library Complexity

#### Read Assignment

Distribution of reads across assigned features, rRNA, unmapped reads, low-quality reads, and filtering categories.

#### Mapping Rate

Percentage of reads assigned to annotated features.

#### Detected Genes

Number of genes/transcripts detected in a sample.

#### Gene Detection Curve

Number of detected transcripts as a function of sequencing depth.

#### Normalized AUC

Area under the gene-detection curve. Higher values indicate more efficient transcript discovery at lower sequencing depth.

#### Tail Gain

Fraction of additional genes detected after ~75% of reads have been sampled. High values suggest sequencing depth has not yet saturated.

#### Top-100 Transcript Fraction

Fraction of reads assigned to the 100 most abundant transcripts. Elevated values may indicate low library complexity.

#### Duplication Rate

Estimated read duplication level.

#### Shannon Entropy

Diversity of transcript abundance distribution.

#### Effective Gene Count

Exponentiated Shannon entropy; interpretable as the effective number of expressed genes.

#### Fragment Length

Estimated library fragment size.

---

## RNA Degradation

#### Gene Body Coverage

Coverage distribution across transcript bodies.

#### Degradation Score

5′/3′ coverage bias used to identify RNA degradation. Values near zero indicate uniform coverage.

---

## Replicate Consistency

#### Sample Correlation

Pearson correlation between samples using highly variable genes.

#### PCA

Principal component analysis of gene expression profiles to identify outliers and batch effects.

#### PCA Distance to Centroid

Within-condition dispersion used to quantify replicate consistency.

#### Top-1000 Gene Jaccard Similarity

Similarity of highly expressed genes between samples.

---

## Differential Expression Diagnostics

#### Size Factors

Normalization factors used for count scaling.

#### Dispersion Estimates

Gene-wise and fitted dispersion estimates used by differential expression methods.

---

## Sequencing QC

#### Per-Base Quality

Base quality scores before and after trimming.

#### Per-Sequence Quality

Distribution of read quality scores.

#### GC Content

GC-content profiles across read positions.

#### Adapter Content

Residual adapter contamination after trimming.

#### Insert Size Distribution

Estimated fragment size distribution for paired-end libraries.
 
