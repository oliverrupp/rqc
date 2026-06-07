# Project Structure & File Format

Detailed guide to setting up your RQC project.

## Directory Hierarchy

```
project_root/
├── organism1/                          # Species/experiment 1
│   ├── reference/
│   │   ├── genome.fa                   # Genome FASTA (REQUIRED)
│   │   ├── annotation.gtf              # Gene annotation (RECOMMENDED)
│   │   ├── annotation.gff3             # OR this (auto-converts to GTF)
│   │   ├── samples.tsv                 # Sample metadata (REQUIRED)
│   │   ├── samples_group1.tsv          # Optional: alternative grouping
│   │   └── [other files created by pipeline]
│   │       ├── transcripts.fa          # Extracted transcripts
│   │       ├── genome.fa.fai           # Genome index
│   │       └── rrna.fa / rRNA.gtf      # rRNA sequences and annotations
│   │
│   ├── reads/                           # Input sequencing data
│   │   ├── sample1_1.fq.gz             # Paired-end forward
│   │   ├── sample1_2.fq.gz             # Paired-end reverse
│   │   ├── sample1_R1.fq.gz            # Alternative naming (R1/R2)
│   │   ├── sample1_R2.fq.gz
│   │   ├── sample2_s.fq.gz             # Single-end reads
│   │   ├── sample3_1.fq.gz
│   │   └── sample3_2.fq.gz
│   │
│   ├── results/                        # Output directory (created by pipeline)
│   │   ├── trimmed/
│   │   ├── salmon/
│   │   ├── salmon_rrna/
│   │   ├── salmon_quantiles/
│   │   ├── falco/
│   │   ├── bam/
│   │   └── index/
│   │
│   ├── report/                         # Generated reports
│   │   ├── samples.report.html
│   │   └── samples_group1.report.html
│   │
│   ├── tsv/
│   │   └── samples/                    # Detailed metrics in TSV format
│   │       ├── complexity.tsv
│   │       ├── pca.tsv
│   │       ├── jaccard.tsv
│   │       └── [many more]
│   │
│   └── benchmark/                      # Execution benchmarks
│       ├── trim_pe.sample1.txt
│       ├── salmon_pe.sample1.txt
│       └── [one per job]
│
├── organism2/
│   └── [same structure]
│
└── .gitignore                          # Exclude large files
```

## File Format Specifications

### Sample Metadata (samples.tsv)

Tab-delimited file with header row. Minimum columns:

```
sample          condition
sample1         treated
sample2         treated
sample3         control
sample4         control
```

**Column specifications:**

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `sample` | string | YES | Unique sample ID. Must match read filename prefix (e.g., `sample1` for `sample1_1.fq.gz`) |
| `condition` | string | YES | Experimental condition for grouping |
| `[batch_factor]` | string | NO | Any additional columns treated as batch factors (e.g., `batch`, `time_point`, `replicate`) |

**Example with batch factors:**

```
sample          condition       batch           time_point
s1_t0_b1        treatment       batch1          T0
s2_t0_b1        treatment       batch1          T0
s3_t6_b1        treatment       batch1          T6
s4_t6_b1        treatment       batch1          T6
s5_t0_b2        control         batch2          T0
s6_t0_b2        control         batch2          T0
```

### Read File Naming

The pipeline detects read files via glob patterns and infers read type from naming:

**Paired-end reads (both patterns supported):**
```
sample_name_1.fq.gz  →  Forward read
sample_name_2.fq.gz  →  Reverse read

OR

sample_name_R1.fq.gz →  Forward read
sample_name_R2.fq.gz →  Reverse read
```

**Single-end reads:**
```
sample_name_s.fq.gz  →  Single-end read
```

**Naming rules:**
- Use consistent naming: either `_1/_2` or `_R1/_R2`, not mixed
- Sample name must not contain `_1`, `_2`, `_s`, `_R1`, `_R2` elsewhere
- File must end with `.fq.gz` (gzip-compressed FASTQ)
- All reads must be gzip-compressed

### Reference Genome

**Format:** FASTA (`.fa`, `.fasta`, or `.fna`)

```
>chr1
ATGCATGCATGCATGC...
>chr2
ATGCATGCATGCATGC...
```

**Requirements:**
- Must be gzip-compressed (`.fa.gz`)
- Sequence names should match annotation file
- Can contain scaffolds/contigs (not limited to chromosomes)

### Gene Annotation

**Supported formats:**

1. **GTF (GFF2) format** - Recommended
   - File: `annotation.gtf`
   - Columns: seqname, source, feature, start, end, score, strand, frame, attributes
   - Must contain `gene_id` and `transcript_id` attributes

2. **GFF3 format** - Automatically converted to GTF
   - File: `annotation.gff3`
   - Will be converted using gffread
   - Parent-child relationships preserved

**Example GTF snippet:**
```
chr1	HAVANA	gene	1000	2000	.	+	.	gene_id "GENE001"; gene_name "ABC1"; gene_type "protein_coding";
chr1	HAVANA	transcript	1000	2000	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; gene_name "ABC1";
chr1	HAVANA	exon	1000	1100	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_number "1";
chr1	HAVANA	exon	1500	2000	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_number "2";
```

### Alternative Sample Groupings

Create multiple sample TSV files for different analyses:

```
reference/
├── samples.tsv              # Primary analysis
├── samples_by_batch.tsv     # Batch-specific analysis
└── samples_by_tissue.tsv    # Tissue-specific analysis
```

The pipeline will generate separate reports for each file:
- `samples.report.html`
- `samples_by_batch.report.html`
- `samples_by_tissue.report.html`

## Directory Validation

The pipeline validates:

✓ **Organism directory exists** and contains:
- `reference/` subdirectory
- `reads/` subdirectory

✓ **Reference directory contains:**
- `genome.fa` file (required)
- Either `annotation.gtf` or `annotation.gff3` (required)
- At least one `samples*.tsv` file (required)

✓ **Reads directory contains:**
- At least one pair of paired-end reads (e.g., `*_1.fq.gz` + `*_2.fq.gz`)
- OR at least one single-end read (e.g., `*_s.fq.gz`)
- All files are `.fq.gz` compressed

✓ **Sample file is valid:**
- Tab-delimited format
- Header row with at least `sample` and `condition` columns
- Sample names match read file prefixes
- No duplicate sample names

Run validation without executing:
```bash
python rqc.py --validate
```

Check which organisms are detected:
```bash
python rqc.py --list-organisms
```

## Best Practices

### Naming Conventions
- Use descriptive, short sample names (avoid special characters)
- Use consistent naming across projects for reproducibility
- Avoid reserved shell characters: `*`, `?`, `[`, `]`, `$`, etc.

### Compression
- Always compress FASTQ files with gzip (`.fq.gz`)
- Use `-9` compression level for long-term storage: `gzip -9 file.fq`
- Smaller files reduce I/O time and storage needs

### Annotation Quality
- Verify genome and annotation are compatible (same version/source)
- Check that sequence names match exactly
- Use primary assemblies when available
- Ensure annotation covers most of the genome

### Sample Design
- Include at least 2-3 biological replicates per condition
- Record and include all relevant batch factors
- Ensure balanced designs where possible
- Document any technical replicates separately

### Storage Organization
```
RNA-seq_projects/
├── project_name_YYYYMM/
│   ├── organism1/
│   │   ├── reference/
│   │   │   ├── genome.fa
│   │   │   ├── annotation.gtf
│   │   │   └── samples.tsv
│   │   └── reads/
│   │       └── [fastq files]
│   └── organisms2_YYYYMM/
└── [future projects]
```

## Example Project Setup

Create a complete minimal example:

```bash
# Create directory structure
mkdir -p myproject/organism1/{reference,reads}

# Download reference (example: small test genome)
wget -O myproject/organism1/reference/genome.fa.gz <url>
wget -O myproject/organism1/reference/annotation.gtf.gz <url>

# Create sample metadata
cat > myproject/organism1/reference/samples.tsv << EOF
sample	condition
S1	control
S2	control
S3	treated
S4	treated
EOF

# Link or copy FASTQ files
ln -s /path/to/reads/S1_1.fq.gz myproject/organism1/reads/
ln -s /path/to/reads/S1_2.fq.gz myproject/organism1/reads/
# ... repeat for all samples

# Run pipeline
python rqc.py myproject/
```

## Troubleshooting

**Error: "No valid organisms found"**
- Check that subdirectories exist: `organism1/reference/`, `organism1/reads/`
- Verify required files: `genome.fa`, annotation file, `samples.tsv`, reads

**Error: "Missing required annotation files"**
- Must have either `.gtf` or `.gff3` file in reference/
- Files should not have additional extensions (e.g., `annotation.gtf.gz` → extract first)

**Error: "No valid samples*.tsv files found"**
- Check samples.tsv file exists in `reference/` directory
- Verify it contains `sample` and `condition` column headers
- Check for tab-delimiters (not spaces)

**Error: "No .fq.gz files found"**
- Ensure reads are in gzip format (`.fq.gz`, not `.fq`)
- Verify files in correct `reads/` subdirectory
- Check naming: `*_1.fq.gz`, `*_2.fq.gz` or `*_R1.fq.gz`, `*_R2.fq.gz` or `*_s.fq.gz`

See [Troubleshooting Guide](TROUBLESHOOTING.md) for more help.
