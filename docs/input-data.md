# Input Data Format

## Directory Layout

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

Multiple `samples*.tsv` are allowd, each file will be analyzed separately. This allows for different projects on the same reference organism.

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
