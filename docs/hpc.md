# Running RQC on HPC Clusters

RQC is designed to run on local workstations as well as HPC environments managed by schedulers such as SLURM, PBS/Torque, LSF, or SGE.

The workflow itself is implemented in Snakemake and can therefore be executed either on a single node or distributed across a cluster.

## Resource Requirements

Resource requirements depend on:

* number of samples
* sequencing depth
* genome size
* selected aligner/quantification tools

Typical requirements for mammalian bulk RNA-seq:

| Step                       | CPU       | Memory  |
| -------------------------- | ----------| ------- |
| QC and trimming            |  8 cores  | 4–16 GB |
| Alignment / quantification | 16 cores  | 8–64 GB |
| Counting                   |  1 cores  | 4–16 GB |
| Report generation          |  1 cores  | 2–8 GB  |

A complete project containing dozens to hundreds of samples is typically executed using cluster scheduling.

---

## Shared Filesystems

RQC assumes that all compute nodes can access:

```text
project/
└── organism1/
      ├── reads/
      ├── reference/
      └── results/
└── organism2/
      └── ...
```

through a shared filesystem such as:

* NFS
* Lustre
* BeeGFS
* GPFS

```bash
rqc --hpc slurm --hpc-config config.yaml
```

---

## Conda Environment

The `rqc` wrapper automatically creates the main required Conda environment for snakemake.

Other Conda environments will be install by snakemake in the `projects` folder.

The first execution may take several minutes while dependencies are installed.

Subsequent executions in the same project reuse existing environments.



---

## Restarting Failed Runs

Snakemake automatically resumes incomplete workflows.

Simply rerun:

```bash
rqc project/ --rerun-incomplete
```

Only missing or outdated outputs will be recomputed.



## Continue pipeline execution on failed jobs

```bash
rqc project/ --keep-going
```
