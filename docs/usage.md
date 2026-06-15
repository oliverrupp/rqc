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
