# Workflow

1. Input validation
2. Read preprocessing and trimming
3. Read alignment or quantification
4. Gene and transcript quantification
5. Library complexity analysis
6. RNA degradation analysis
7. Replicate consistency assessment
8. Differential expression diagnostics
9. Report generation

## Output

For every sample table (`samples*.tsv`), the pipeline generates:

```text
organism/
└── report/
    └── samples.report.html
```

Interactive reports contain plots, tables, PCA visualizations, QC summaries, and sequencing diagnostics.
