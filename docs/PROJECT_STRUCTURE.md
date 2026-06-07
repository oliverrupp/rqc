# Project Structure

## Directory Layout

The RQC pipeline expects your project to be organized as follows:

```
project_dir/
├── rqc.py                           # This wrapper script
├── rqc.smk                          # Snakemake workflow
├── subproject1/
│   ├── reference/
│   │   ├── genome.fa                # Reference genome (FASTA)
│   │   ├── annotation.gtf           # Gene annotations (GTF)
│   │   └── samples.tsv              # Sample metadata
│   └── reads/
│       ├── sample1_1.fq.gz
│       ├── sample1_2.fq.gz
│       ├── sample2_1.fq.gz
│       ├── sample2_2.fq.gz
│       └── ...
├── subproject2/
│   ├── reference/
│   │   ├── genome.fa
│   │   ├── annotation.gtf
│   │   └── samples.tsv
│   └── reads/
│       └── ...
└── subproject3/
    └── ...
```

## Required Files per Subproject

### 1. Reference Files

#### `reference/genome.fa`
- **Format**: FASTA (nucleotide sequences)
- **Content**: Reference genome sequences
- **Example**:
  ```fasta
  >chr1
  ATCGATCGATCG...
  >chr2
  GATCGATCGATC...
  ```

#### `reference/annotation.gtf`
- **Format**: GTF (Gene Transfer Format)
- **Content**: Gene/transcript annotations
- **Example**:
  ```
  chr1	ensembl	gene	1000	2000	.	+	.	gene_id "ENSG00000001"; gene_name "GENE1"
  chr1	ensembl	transcript	1000	2000	.	+	.	gene_id "ENSG00000001"; transcript_id "ENST0000001"
  ```

#### `reference/samples.tsv`
- **Format**: Tab-separated values (TSV)
- **Required Columns**: `sample`, `condition` (at minimum)
- **Optional Columns**: batch, replicate, or other metadata
- **Example**:
  ```
  sample	condition	batch	replicate
  sample1	control	batch1	1
  sample2	control	batch1	2
  sample3	treated	batch1	1
  sample4	treated	batch1	2
  ```

### 2. Read Files

#### `reads/` Directory
- **Format**: Gzipped FASTQ files (`.fq.gz` or `.fastq.gz`)
- **Naming Convention**: `{sample_name}_{read_type}.fq.gz`
- **Read Types**:
  - `_1.fq.gz` - Forward reads (paired-end)
  - `_2.fq.gz` - Reverse reads (paired-end)
  - `_s.fq.gz` - Single-end reads
  - `_R1.fq.gz`, `_R2.fq.gz` - Alternative paired-end naming
  - `_S.fq.gz` - Alternative single-end naming

**Example File Structure**:
```
reads/
├── sample1_1.fq.gz
├── sample1_2.fq.gz
├── sample2_1.fq.gz
├── sample2_2.fq.gz
├── sample3_1.fq.gz
├── sample3_2.fq.gz
└── sample4_1.fq.gz
```

## Validation Rules

The RQC pipeline automatically validates:

1. **Project directory exists** and is readable
2. **Each subproject contains**:
   - `reference/genome.fa` file
   - `reference/annotation.gtf` file
   - `reference/samples.tsv` file with `sample` and `condition` columns
   - `reads/` directory with at least one `.fq.gz` file
3. **Sample names** in `samples.tsv` match sequencing file names

## Multiple Subprojects

You can organize multiple independent analyses in the same project:

```bash
# Run only subproject1
rqc.py --subproject subproject1 --max-cpus 16

# Run multiple subprojects
rqc.py --subproject subproject1,subproject2 --max-cpus 16

# Run all subprojects (default)
rqc.py --max-cpus 16

# List available subprojects
rqc.py --list-subprojects
```

Each subproject will generate its own report and output files.

## Best Practices

1. **Use descriptive names** for samples and subprojects
2. **Keep reference files up-to-date** and consistent across runs
3. **Include all relevant metadata** in `samples.tsv` for better QC interpretation
4. **Use consistent naming** for read files (e.g., all paired-end with `_1`/`_2` suffix)
5. **Verify file integrity** before running (check gzip files are not corrupted)