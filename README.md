This repository is a rewritten fork of the official PacBio HiFi-16S workflow, redesigned for improved performance, modularity, and flexibility on HPC systems.

# Overview

This pipeline processes PacBio HiFi 16S amplicon sequencing data using a modular Nextflow DSL2 workflow. It includes quality control, filtering, denoising with DADA2, and taxonomic assignment.

The refactor focuses on:

- Improved parallelisation

- Separation of pipeline stages

- Reproducibility and benchmarking

- Compatibility with modern PacBio Revio data

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


# Pipeline Structure

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

# Usage

## Test pipeline

```
nextflow run main.nf -profile test,docker # or signularity, or conda
```

## Basic Command

```
nextflow run main.nf \
  --input samplesheet.tsv \
  --meta_data meta_data.tsv
  --outdir results \
  -profile conda

```
