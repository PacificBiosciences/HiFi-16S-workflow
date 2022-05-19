# pb-16S-analysis Pipeline using QIIME 2

## Installation and usage
This pipeline runs using Nextflow (Version 22 and above). All softwares dependencies are
managed via `Conda`. We recommend installing [`mamba`](https://github.com/mamba-org/mamba)
to speed up the conda environment installation. The default `nextflow.config` file 
enables the use of `mamba` by default. You can install Nextflow following the instruction
from Nextflow [documentation](https://www.nextflow.io/docs/latest/getstarted.html) or via Conda:

```
conda install -c bioconda nextflow
```

The taxonomy classification step of the pipeline requires a few databases. These can be downloaded
from the following link:

- `--vsearch_db` and `--vsearch_tax`
  - [`silva-138-99-seqs.qza`](https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza)
  - [`silva-138-99-tax.qza`](https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza)
- `--silva_db`
  - [`silva_nr99_v138.1_wSpecies_train_set.fa.gz`](https://zenodo.org/record/4587955)
- `--gtdb_db`
  - [`GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz`](https://zenodo.org/record/4735821)
- `--refseq_db`
  - [`RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz`](https://zenodo.org/record/4735821)

After installing Nextflow and `mamba`(Optional but recommended), clone the repository:

```
git clone https://github.com/proteinosome/pb-16S-nf.git
```

This repository contains a set of small test data as well as example samples and metadata
TSV file that can be used to test the pipeline. Note that you need to change the
path of the FASTQ in `test_samples.tsv` to point to the correct location of the
test dataset.

After cloning the repository, run the following command in the cloned folder
to see the options for the pipeline:

```
nextflow run main.nf --help

Usage:
  This pipeline takes in the standard sample manifest and metadata file used in
  QIIME 2 and produces QC summary, taxonomy classification results and visualization.

  For samples TSV, two columns named "sample-id" and "absolute-filepath" are
  required. For metadata TSV file, at least two columns named "sample_name" and
  "condition" to separate samples into different groups.

  nextflow run main.nf --input samples.tsv --metadata metadata.tsv \
    --dada2_cpu 8 --vsearch_cpu 8

  By default, sequences are first trimmed with cutadapt (higher rate compared to using DADA2
  ) using the example command below. You can skip this by specifying "--skip_primer_trim"
  if the sequences are already trimmed. The primer sequences used are the F27 and R1492
  primers for full length 16S sequencing.

  Other important options:
  --filterQ    Filter input reads above this Q value (default: 30).
  --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
              this value will be discarded (default: 2)
  --min_len    Minimum length of sequences to keep (default: 1000)
  --max_len    Maximum length of sequences to keep (default: 1600)
  --pooling_method    QIIME 2 pooling method for DADA2 denoise (default: "pseudo"),
                      see QIIME 2 documentation for more details
  --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 5)
  --vsearch_identity    Minimum identity to be considered as hit (default 0.97)
  --rarefaction_depth    Rarefaction curve "max-depth" parameter. By default the pipeline
                         automatically select a cut-off above the minimum of the denoised
                         reads for >90% of the samples. This cut-off is stored in a file called
                         "rarefaction_depth_suggested.txt" file in the results folder
                         (default: null)
  --dada2_cpu    Number of threads for DADA2 denoising (default: 8)
  --vsearch_cpu    Number of threads for VSEARCH taxonomy classification (default: 8)
  --cutadapt_cpu    Number of threads for primer removal using cutadapt (default: 16)
  --outdir    Output directory name (default: "results")
  --vsearch_db  Location of VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Location of VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --silva_db   Location of Silva 138 database for taxonomy classification
  --gtdb_db    Location of GTDB r202 for taxonomy classification
  --refseq_db    Location of RefSeq+RDP database for taxonomy classification
  --skip_primer_trim    Skip all primers trimming (switch off cutadapt and DADA2 primers
                        removal) (default: trim with cutadapt)
```

To test the pipeline, run this example below. Note that the path of the database needs
to be changed to the location on your server. Conda environment will by default be created at 
`$HOME/nf_conda` folder unless changed in the `nextflow.config` file. Once the conda environment
is created it will be reused by any future run.

```
nextflow run main.nf --input test_sample.tsv \
    --metadata test_metadata.tsv -profile conda \
    --dada2_cpu 32 --vsearch_cpu 32 --outdir results \
    --vsearch_db /path/to/silva-138-99-seqs.qza \
    --vsearch_tax /path/to/silva-138-99-tax.qza \
    --silva_db /path/to/silvaDB \
    --gtdb_db /path/to/gtdbDB \
    --refseq_db /path/to/refseqDB
```

To run this pipeline on your data, create the sample TSV and metadata TSV following
the test data format (For metadata, if you do not have any grouping, you can just
put any words in the "condition" column) and run the workflow similar to the above. 
Remember to specify the `--outdir` directory to avoid overwriting existing results.

The pipeline uses Slurm scheduler by default to run jobs on HPC. This can be changed
in the `nextflow.config` file under `executor`. See Nextflow 
[documentation](https://www.nextflow.io/docs/latest/executor.html) on the available
executors. CPUs for `VSEARCH`, `DADA2` and `cutadapt` can be specified as command line
parameters as shown above. For all the other processes, they use any of the default
labels in `nextflow.config` and can be changed according to your need.

Pipeline is still under active development. The nextflow.config file by default will 
generate workflow DAG and resources report to help benchmarking the resources
required.

## Outputs
In the output folder, you will find many useful results. In the `results` folder,
there is a HTML file named `visualize_biome.html` that provides an overview report
of the 16S community with the taxonomy tables generated for filtering in web browser.
The HTML file in `krona_html` can be opened directly in the browser for Krona plot
visualization of the community. All ASV sequences in FASTA format can be found
in the folder `dada2` as `dada2_ASV.fasta`.

All the files ending with `.qzv` extension can be
opened in [QIIME 2 View](https://view.qiime2.org/) directly for interactive 
visualization, e.g. `taxa_barplot.qzv` for taxonomy barplot. For advanced users,
`merged_freq_tax.tsv` contains the OTU count matrix that can be loaded into any
language of choice for further processing or plots. You may also import the
BIOM format file `feature-table-tax.biom` with many packages such as the popular
`phyloseq`.

## Frequently asked questions (FAQ)
* A lot of my reads are lost in the denoise stage, what's going on?

This can happen in extremely diverse community such as soil where the ASVs are of very low abundance.
In each sample, the reads supporting the ASV are very low and may not pass DADA2 threshold to qualify
as a cluster. See [here](https://github.com/benjjneb/dada2/issues/841) for a discussion on
the algorithm.

* I'm getting `Conda` "Safety" error indicating corrupted package or that some
pipeline steps are not able to find specific command line tools (e.g. qiime).

Sometimes the conda cache can become corrupted if you run many workflows
in parallel before the environment was created, thus causing conflicts between
the different workflow competing to create the same environment. We recommend running
the test dataset above and wait for it to finish first so the conda
environment is created successfully. Subsequent runs will use the same
environment and will not need to recreate them. Nevertheless, if the errors
happen, try running `conda clean -a` and remove the offending `conda` packages cache
in the cache directory (e.g. if the erorr happens for QIIME, delete any folder in `conda info` 
"package cache" containing QIIME).

You can try to install the QIIME 2 environment directly to inspect any error
messages:

`mamba env create -n q2_test -f qiime2-2022.2-py38-linux-conda.yml`

* I've received/downloaded 16S FASTQ that already has the primers trimmed, can I skip
primers removal?

We recommend using the pipeline to trim the primers as it works well for HiFi sequencing
data. However, there are many public dataset that may already have the full length
primers trimmed, in which case you can specify `--skip_primer_trim` to skip
primers trimming. If unsure, run with default pipeline and the cutadapt demultiplexing rate
(in the file `results/samples_demux_rate.tsv`) should be close to zero for all
samples if the primers are already trimmed.

* How does the taxonomy classification work?

We use the `assignTaxonomy` function from DADA2 to carry out taxonomy assignment.
This pipeline uses 3 databases to classify the ASVs and the priority of assignment
is Silva 138, followed by GTDB r202, then lastly RefSeq + RDP. This means for example
if an ASV is not assigned at Species level using Silva, it will check if it can be assigned
with GTDB. This ensure we assign as many ASVs as possible. There is also a VSEARCH
taxonomy classification using Silva database only in the file called `results/vsearch_merged_freq_tax.tsv`
that may work well in some cases as an alternative.

## DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, 
ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES 
OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY 
WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS 
FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE 
OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR 
APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY 
OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES 
DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
