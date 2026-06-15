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

Align reads to genome with STAR, use the aligned reads for quantification

```bash
rqc --alignment
```



