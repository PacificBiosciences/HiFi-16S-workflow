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

# Reference databases

This workflow supports taxonomic classification of PacBio HiFi amplicon
sequence variants (ASVs) using reference databases appropriate for the
sequenced marker.

A deliberate design principle of the workflow is that **only officially
published reference databases or officially published classifier-specific
derivatives are used**. The workflow does not generate custom SINTAX
databases, reconstruct taxonomies from separate sequence and taxonomy files,
or convert databases between classifier formats.

This keeps the taxonomic references reproducible and avoids introducing
pipeline-specific decisions into the construction of reference databases.

## Supported databases

The appropriate reference database depends on the marker being sequenced.

| Database       | Marker         | Typical target                 | Naive Bayes | VSEARCH-LCA |
| -------------- | -------------- | ------------------------------ | ----------: | ----------: |
| SILVA          | 16S rRNA       | Bacteria and Archaea           |           ✓ |           — |
| GTDB           | 16S rRNA       | Bacteria and Archaea           |           ✓ |           — |
| Greengenes2    | 16S rRNA       | Bacteria and Archaea           |           ✓ |           — |
| EUKARYOME SSU  | 18S / SSU rRNA | Eukaryotes                     |           ✓ |           ✓ |
| EUKARYOME LSU  | 28S / LSU rRNA | Eukaryotes                     |           ✓ |           ✓ |
| EUKARYOME ITS  | ITS            | Eukaryotes, particularly fungi |           ✓ |           ✓ |
| EUKARYOME long | SSU–ITS–LSU    | Long eukaryotic amplicons      |           ✓ |           ✓ |

The database identifiers used by the workflow are:

```text
silva
gtdb
gg2
euk_ssu
euk_lsu
euk_its
euk_long
```

## 16S rRNA data

For bacterial and archaeal 16S rRNA amplicons, the workflow supports three
reference databases:

* **SILVA**
* **GTDB**
* **Greengenes2**

These databases are currently used with the Naive Bayes taxonomy
classification workflow.

### SILVA

SILVA provides a curated collection of aligned small- and large-subunit
ribosomal RNA sequences across all domains of life. In this workflow, the
SILVA SSU reference is used for taxonomic classification of bacterial and
archaeal 16S rRNA sequences.

The currently configured classifier is based on **SILVA 138.2**.

### GTDB

The Genome Taxonomy Database (GTDB) provides a standardized genome-based
taxonomy for Bacteria and Archaea.

The workflow uses an officially published DADA2-compatible SSU reference
derived from **GTDB release R220**. GTDB is particularly useful when a
genome-based bacterial and archaeal taxonomy is desired.

### Greengenes2

Greengenes2 provides an updated reference taxonomy integrating microbial
genome and marker-gene information.

The workflow currently uses the officially published
**Greengenes2 2024.09** DADA2-compatible reference.

### Why VSEARCH-LCA is not used for these databases

The workflow does not construct its own VSEARCH/SINTAX versions of SILVA,
GTDB, or Greengenes2.

Although it is technically possible to convert sequence and taxonomy files
into a SINTAX-compatible FASTA, doing so would require the workflow to define
its own rules for taxonomy parsing, rank normalization, identifier matching,
and reference construction.

Instead, the workflow uses the officially published classifier-ready
references directly. Consequently, SILVA, GTDB, and Greengenes2 are currently
available through the Naive Bayes classification branch only.

## Eukaryotic marker data

Eukaryotic amplicons are classified using the **EUKARYOME** reference
database.

EUKARYOME provides marker-specific reference sets for several commonly used
eukaryotic amplicons. Importantly, EUKARYOME publishes databases specifically
prepared for different taxonomic classification algorithms.

The workflow therefore uses the official EUKARYOME DADA2 references for the
Naive Bayes classifier and the official EUKARYOME SINTAX references for
VSEARCH.

No conversion between these formats is performed by the workflow.

### 18S / SSU

For 18S rRNA amplicon sequencing, use:

```text
euk_ssu
```

This corresponds to the EUKARYOME small-subunit (SSU) reference database.

Typical input consists of eukaryotic 18S rRNA amplicons. Both classification
methods are supported:

```text
euk_ssu/
├── nb/
│   └── DADA2_EUK_SSU_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_SSU_v2.1.fa.gz
```

### 28S / LSU

For 28S rRNA amplicon sequencing, use:

```text
euk_lsu
```

This corresponds to the EUKARYOME large-subunit (LSU) reference database.

Both Naive Bayes and VSEARCH-LCA classification are supported:

```text
euk_lsu/
├── nb/
│   └── DADA2_EUK_LSU_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_LSU_v2.1.fa.gz
```

### ITS

For eukaryotic ITS amplicons, use:

```text
euk_its
```

ITS is widely used for eukaryotic taxonomic profiling and is particularly
important for fungal identification.

Both official EUKARYOME classifier formats are supported:

```text
euk_its/
├── nb/
│   └── DADA2_EUK_ITS_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_ITS_v2.1.fa.gz
```

### Long SSU–ITS–LSU amplicons

PacBio HiFi sequencing makes it possible to sequence long eukaryotic
amplicons spanning multiple ribosomal markers. For combined
SSU–ITS–LSU sequences, use:

```text
euk_long
```

This reference should be used when the sequenced amplicon spans the extended
eukaryotic ribosomal region rather than an isolated SSU, LSU, or ITS marker.

Both classifier formats are available:

```text
euk_long/
├── nb/
│   └── DADA2_EUK_long_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_long_v2.1.fa.gz
```

The long reference should not be selected simply because the sequencing
technology is PacBio HiFi. Database selection is determined by the
**amplified marker region**, not by read length or sequencing platform.

For example, a full-length 18S amplicon generated with PacBio HiFi should
still use `euk_ssu`, whereas an amplicon designed to span SSU, ITS, and LSU
should use `euk_long`.

## Database selection

In general, the marker determines the appropriate database:

```text
16S
 ├── silva
 ├── gtdb
 └── gg2

18S / SSU
 └── euk_ssu

28S / LSU
 └── euk_lsu

ITS
 └── euk_its

SSU–ITS–LSU
 └── euk_long
```

Databases targeting different marker regions should not be treated as
interchangeable.

## Classifier-specific references

The workflow keeps classifier-specific reference files in separate
directories:

```text
<database>/
├── nb/
└── vsearch/
```

`nb/` contains references intended for the Naive Bayes/DADA2 taxonomy
classification workflow.

`vsearch/` contains official SINTAX-formatted references suitable for
VSEARCH-based taxonomy assignment and LCA classification.

For example:

```text
euk_ssu/
├── nb/
│   └── DADA2_EUK_SSU_v2.1.fa.gz
└── vsearch/
    └── SINTAX_EUK_SSU_v2.1.fa.gz
```

This separation is intentional. The two files are independently published
classifier-specific representations of the reference database and should not
be generated from one another by the pipeline.

## Database download

All supported databases can be installed through the workflow:

```bash
nextflow run main.nf \
    --download_db true \
    --download_targets silva,gtdb,gg2,euk_ssu,euk_lsu,euk_long,euk_its
```

The location of the installed databases is controlled by `--db_base_dir`.

The database manifest in:

```text
conf/databases.yml
```

defines the database versions, filenames, marker types, and official download
locations.

For example, an installation may have the following structure:

```text
databases/
├── silva/
│   └── nb/
├── gtdb/
│   └── nb/
├── gg2/
│   └── nb/
├── euk_ssu/
│   ├── nb/
│   └── vsearch/
├── euk_lsu/
│   ├── nb/
│   └── vsearch/
├── euk_its/
│   ├── nb/
│   └── vsearch/
└── euk_long/
    ├── nb/
    └── vsearch/
```

Downloaded archives are treated only as installation artifacts. The final
database directories contain the decompressed and prepared reference FASTA
files required by the corresponding classifier.

## Reference database policy

Reference database preparation can have a substantial effect on taxonomic
classification. Decisions such as removing sequences, translating taxonomic
ranks, resolving missing ranks, modifying identifiers, or constructing
classifier-specific headers can change classification results.

For this reason, this workflow follows a conservative database policy:

> **Classifier-specific reference databases are obtained from officially
> published resources whenever available and are not reconstructed by the
> workflow.**

In particular:

* SILVA, GTDB, and Greengenes2 use published DADA2-compatible references.
* EUKARYOME uses the officially published DADA2 references for Naive Bayes
  classification.
* EUKARYOME uses the officially published SINTAX references for VSEARCH
  classification.
* The workflow does not generate SINTAX databases from SILVA, GTDB, or
  Greengenes2.
* The workflow does not maintain a custom taxonomy-to-SINTAX conversion
  procedure.
* Download and extraction steps do not alter the biological content or
  taxonomy of the published reference database.

This design makes the provenance of each taxonomic reference explicit and
keeps database construction separate from the analysis workflow.

