This repository is a rewritten fork of the official PacBio HiFi-16S workflow, redesigned for improved performance, modularity, and flexibility on HPC systems.

# Overview

This pipeline processes PacBio HiFi 16S amplicon sequencing data using a modular Nextflow DSL2 workflow. It includes quality control, filtering, denoising with DADA2, and taxonomic assignment. The main output of the pipeline is a count table with taxonomic assignment.

Please create an issue for bugs and feature requests.

The refactor focuses on:

- Improved parallelisation

- Separation of pipeline stages

- Reproducibility and benchmarking

- Compatibility with modern PacBio Revio data

- Database management


# Key Features

**Parallelised Execution**

- Parallelisation handled at the Nextflow level, not inside R

- Improved scalability on HPC systems

**Modular DADA2 Workflow**

- Each step is split into independent processes:

  - Filtering

  - Error learning

  - Denoising

  - Sequence table construction

  - Chimera removal

Easier debugging and benchmarking

**Independent Sample Processing**

- Trimming, filtering, and denoising operate on individual FASTQ files

- Enables:

  - Fine-grained parallelisation

  - Failure isolation per sample


**Independent Denoising Mode**

- Each library can be denoised separately (independent mode)

- Useful for:

  - Benchmarking
  - Large-scale datasets

**Note**: May slightly affect ASV consistency across samples

**Improved Error Model Handling**

- Supports binned quality scores (PacBio Revio)

- Error models can be:

  - Learned once

  - Reused across runs (⚠️ use with caution—see below)

**Flexible Taxonomic Assignment**

- Taxonomy assignment split into one process per database

- Supports multiple databases (e.g. SILVA, GTDB, Greengenes2)

- Enables:

  - Easy comparison across databases
  
  - Independent execution

# Missing Features

- Denoising with pooling strategy

- Phyloseq creation as final output

- Small Report generation

# Pipeline Structure

```
Input FASTQ
   ↓
QC + Filtering
   ↓
Primer trimming
   ↓
Error model learning
   ↓
Denoising (independent)
   ↓
Sequence table construction
   ↓
Chimera removal
   ↓
ASV filtering
   ↓
Taxonomic assignment (per database)
```

# Usage

## Test pipeline

```
nextflow run main.nf -profile test,docker # or signularity, or conda
```
## Download Databases

Make sure to set the preferred location of your databases in the `nextflow.config`.

```
nextflow run main.nf \
  --download_db \
  --download_targets silva,gg2,gtdb \
  -profile conda
```

The structure of the database directory looks like this:

```
├── gg2
│   ├── nb
│   │   └── gg2_2024_09_toSpecies_trainset.fa.gz
│   └── vsearch
│       ├── sequences.fasta
│       └── taxonomy.tsv
├── gtdb
│   ├── nb
│   │   └── GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz
│   └── vsearch
│       ├── sequences.fasta
│       └── taxonomy.tsv
└── silva
    ├── nb
    │   └── silva_nr99_v138.2_toSpecies_trainset.fa.gz
    └── vsearch
        ├── sequences.fasta
        └── taxonomy.tsv
```

## Basic Command

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --meta_data meta_data.tsv
  --outdir results \
  -profile conda

```

## Cluster execution

You can execute the pipeline on a cluster by using either a custom nextflow profile, or 
by adding a custom config. See `custom_slurm.config` for an example.

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --metadata meta_data.tsv
  --outdir results \
  -profile conda
  -c custom_slurm.config
```


