# Workflow

1. Input validation
2. Read preprocessing and trimming
3. Optional: Read alignment 
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
      └── samples*.report.html
└── tsv/samples*/
          ├── g_TPM.tsv
	      ├── t_TPM.tsv
          ├── g_TMM.tsv
	      ├── t_TMM.tsv
		  └── ...
```

Interactive reports contain plots, tables, PCA visualizations, QC summaries, and sequencing diagnostics.
Raw and noramlized counts will be available in the `tsv/samples*` folders. 

### Available Counts
All counts are available on gene-level (`g_`) and transcript-level (`t_`).

- raw counts (estimated values by salmon)
- TPM
- TMM
- geTMM

All normalizations will be computed separately for each `samples*.tsv` file.
