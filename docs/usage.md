# Usage

The main entry point is the `rqc` shell script, which automatically creates and manages the required Conda environment and executes the pipeline.

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

Run a genome-guided assembly of the RNA-seq reads, if no annotation is available.

```bash
rqc --assembly
```

Run [BUSCO](https://busco.ezlab.org/) analysis with specific lineage

```bash
rqc --busco LINEAGE
```


Snakemake automatically resumes incomplete workflows.
Only missing or outdated outputs will be recomputed.

Simply rerun:

```bash
rqc --rerun-incomplete
```


Continue pipeline execution on failed jobs

```bash
rqc --keep-going
```
