# HiFi Full-length 16S analysis with pb-16S-nf

- Table of Contents
  - [Workflow overview and output](#workflow-overview-and-output)
  - [Installation and usage](#installation-and-usage)
  - [HPC and job scheduler usage](#hpc)
  - [Speeding up the denoising process](#pooling)
  - [Run time and compute requirements](#runtime)
  - [Frequently asked questions (FAQ)](#faq)
  - [References](#references)
  - [DISCLAIMER](#disclaimer)

## The pipeline is currently under active development; we welcome your feedback to help improve it

## Workflow overview and output

![alt text](misc/pipeline_workflow.png "Workflow diagram")

This Nextflow pipeline is designed to process PacBio HiFi full-length 16S data into high- quality amplicon sequence variants (ASVs) using `QIIME 2` and `DADA2`.
The pipeline provides a set of visualizations using the `QIIME 2` framework for interactive plotting. The pipeline generates an HTML report for
the important statistics and top taxonomies. The outputs and stages of this pipeline are documented [here](pipeline_overview.md).

We provide a sample report generated using this pipeline based on 8 replicates from the ATCC MSA-1003
mock community sequenced on a Sequel II system ([Link](https://www.pacb.com/connect/datasets/)). Right-click this
[link](examples/results/visualize_biom.html?raw=1), save it on
your computer, then double-click to open the sample report. All other important outputs from the pipeline are available
in the [`examples`](examples) folder when you clone this repository.

## Installation and usage

This pipeline runs using Nextflow Version 22 and later. If you have Singularity or Docker on your
cluster, we recommend using Singularity or Docker to run the pipeline by specifying `-profile singularity` or
`-profile docker` when running the pipeline. Singularity will pull the docker images to the folder `$HOME/nf_conda/singularity`.

By default, all software dependencies are managed using `Conda`. Nextflow will use `Conda` to build
the required environment so there is no need to manually build environments.
You can install Nextflow by following the steps here ([documentation](https://www.nextflow.io/docs/latest/getstarted.html))
or by using `Conda` itself:

```
conda install -c bioconda nextflow

# If this is your first time using conda
conda init
```

After installing Nextflow, clone the repository and download databases using the following commands. To update the pipeline in the future,
type `git pull`.

```
git clone https://github.com/PacificBiosciences/HiFi-16S-workflow.git
cd HiFi-16S-workflow
nextflow run main.nf --download_db
# With singularity (If you use singularity, add -profile singularity to all Nextflow-related command)
nextflow run main.nf --download_db -profile singularity
```

After downloading the databases, run the following command in the cloned folder to see the options for the pipeline:

```
nextflow run main.nf --help

  Usage:
  This pipeline takes in the standard sample manifest and metadata file used in
  QIIME 2 and produces QC summary, taxonomy classification results and visualization.

  For samples TSV, two columns named "sample-id" and "absolute-filepath" are
  required. For metadata TSV file, at least two columns named "sample_name" and
  "condition" to separate samples into different groups.

  nextflow run main.nf --input samples.tsv --metadata metadata.tsv \\
    --dada2_cpu 8 --vsearch_cpu 8

  By default, sequences are first trimmed with cutadapt. If adapters are already trimmed, you can skip 
  cutadapt by specifying "--skip_primer_trim".

  Other important options:
  --front_p    Forward primer sequence. Default to F27. (default: AGRGTTYGATYMTGGCTCAG)
  --adapter_p    Reverse primer sequence. Default to R1492. (default: AAGTCGTAACAAGGTARCY)
  --filterQ    Filter input reads above this Q value (default: 20).
  --downsample    Limit reads to a maximum of N reads if there are more than N reads (default: off)
  --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
              this value will be discarded (default: 2)
  --minQ    DADA2 minQ parameter. Reads with any base lower than this score 
            will be removed (default: 0)
  --min_len    Minimum length of sequences to keep (default: 1000)
  --max_len    Maximum length of sequences to keep (default: 1600)
  --pooling_method    QIIME 2 pooling method for DADA2 denoise see QIIME 2 
                      documentation for more details (default: "pseudo", alternative: "independent") 
  --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --min_asv_totalfreq    Total frequency of any ASV must be above this threshold
                         across all samples to be retained. Set this to 0 to disable filtering
                         (default 5)
  --min_asv_sample    ASV must exist in at least min_asv_sample to be retained. 
                      Set this to 0 to disable. (default 1)
  --vsearch_identity    Minimum identity to be considered as hit (default 0.97)
  --rarefaction_depth    Rarefaction curve "max-depth" parameter. By default the pipeline
                         automatically select a cut-off above the minimum of the denoised 
                         reads for >80% of the samples. This cut-off is stored in a file called
                         "rarefaction_depth_suggested.txt" file in the results folder
                         (default: null)
  --dada2_cpu    Number of threads for DADA2 denoising (default: 8)
  --vsearch_cpu    Number of threads for VSEARCH taxonomy classification (default: 8)
  --cutadapt_cpu    Number of threads for primer removal using cutadapt (default: 16)
  --outdir    Output directory name (default: "results")
  --vsearch_db Location of VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Location of VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --silva_db   Location of Silva 138 database for taxonomy classification 
  --gtdb_db    Location of GTDB r202 for taxonomy classification
  --gg2_db    Location of GreenGenes2 database for taxonomy classification
  --db_to_prioritize    Prioritize this database for taxonomy assignment. Must be one of the 
                        following: GG2, GTDB, Silva (default: GG2)
  --skip_primer_trim    Skip all primers trimming (switch off cutadapt and DADA2 primers
                        removal) (default: trim with cutadapt)
  --skip_nb    Skip Naive-Bayes classification (only uses VSEARCH) (default: false)
  --colorby    Columns in metadata TSV file to use for coloring the MDS plot
               in HTML report (default: condition)
  --run_picrust2    Run PICRUSt2 pipeline. Note that pathway inference with 16S using PICRUSt2
                    has not been tested systematically (default: false)
  --download_db    Download databases needed for taxonomy classification only. Will not
                   run the pipeline. Databases will be downloaded to a folder "databases"
                   in the Nextflow pipeline directory.
  --publish_dir_mode    Outputs mode based on Nextflow "publishDir" directive. Specify "copy"
                        if requires hard copies. (default: symlink)
  --version    Output version
```

To test the pipeline, run the example below. Note that the database paths should be changed to their respective locations on your server if they are different. (See the parameters above.) If you
follow the command above, the databases will be downloaded into a `databases` folder in the `pb-16S-nf` folder
and you do not need to specify the path. The conda environment will be created by default in the
`$HOME/nf_conda` folder unless changed in the `nextflow.config` file. Once the conda environment is created, it will be reused by any future run.

```
# Create sample TSV for testing
echo -e "sample-id\tabsolute-filepath\ntest_data\t$(readlink -f test_data/test_1000_reads.fastq.gz)" > test_data/test_sample.tsv

nextflow run main.nf --input test_data/test_sample.tsv \
    --metadata test_data/test_metadata.tsv -profile conda \
    --outdir results

# To test using Singularity or docker (change singularity to docker)
nextflow run main.nf --input test_data/test_sample.tsv \
    --metadata test_data/test_metadata.tsv -profile singularity \
    --outdir results
```

To run this pipeline on your data, create the sample TSV and metadata TSV following
the test data format (for metadata, if you do not have any grouping, put any words in the "condition" column)
and run the workflow similar to the above.
Remember to specify the `--outdir` directory to avoid overwriting existing results.

## HPC and job scheduler usage <a name="hpc"></a>

The pipeline uses "Local" by default to run jobs on HPC. This can be changed
in the `nextflow.config` file under `executor` to use HPC schedulers such as
Slurm, SGE and so on using Nextflow's native support. For example, to use Slurm,
open the `nextflow.config` file and change `executor = 'Local'` to `executor = 'slurm'` and
specify the partition to be used using `queue = PARTITION`. See the Nextflow
[documentation](https://www.nextflow.io/docs/latest/executor.html) for the available
executors and parameters. CPUs for `VSEARCH`, `DADA2` and `cutadapt` can be specified as command-line
parameters. For all the other processes, they use any of the default
labels in `nextflow.config` and can be changed according to your needs.

Note that the `nextflow.config` file, by default, will
generate the workflow DAG and resources reports to help in benchmarking the resources
required. See the `report_results` folder created after the pipeline finishes running
for the DAG and resources report.

## Speeding up `DADA2` denoise <a name="pooling"></a>

By default, the pipeline pools all samples into one single `qza` file for `DADA2` denoise (using
the default pseudo-pooling approach by `DADA2`). This is designed to maximize the sensitivity to
low frequency ASVs. For example, an ASV with just 2 reads in sample 1 may be discarded, but if the same
exact ASV is seen in another sample, this gives the algorithm higher confidence that it is real.
However, when the samples are highly diverse (such as with environmental samples), this can become very slow.

If a (possibly) minor loss in sensitivity is acceptable, the pipeline allows you to "split" the
input samples into different groups that will be denoised separately. This is done using a `pool`
column in the `metadata.tsv` input. Example:

```
sample_name     condition       pool
bc1005_bc1056   RepA    RepA
bc1005_bc1057   RepA    RepA
bc1005_bc1062   RepA    RepA
bc1005_bc1075   RepA    RepA
bc1005_bc1100   RepB    RepB
bc1007_bc1075   RepB    RepB
bc1020_bc1059   RepB    RepB
bc1024_bc1111   RepB    RepB
```

The TSV above will split the 8 samples into two groups (RepA and RepB) and denoise them separately.
After denoising, all denoised ASVs and statistics are merged again for downstream filtering and
processing. This allows you to maximize sensitivity *within* a group of samples and speed up
the pipeline considerably. On the other hand, if each sample has been sequenced deeply, you can
denoise each sample *individually* by setting a unique group for each sample (e.g. replicating
the `sample_name` column as the `pool` column) to process the samples quickly.

Note that this is different from the `--pooling_method` parameter. The `--pooling_method` parameter
controls how the samples are pooled *within* the `DADA2` denoise step. The `pool` column in the metadata
TSV file on the other hand splits each sample into separate `qza` files according to the pool annotation. E.g.
if you have 8 samples and set the `pool` column into two groups, you will have 2
separate `qza` files for the `DADA2` denoise step. Within each of the two denoising steps, the `--pooling_method`
parameter controls whether DADA2 splits the samples for denoising or not.

## Run time and compute requirements <a name="runtime"></a>

We recommend at least 32 CPUs for most sample types. The run time highly depends on the
complexity of the samples in addition to the total number of reads. Shown here are examples of
run times for data tested with this pipeline using 32 CPUs:

|     Sample types      |     Number of samples    |     Number of FL reads    |     Total ASVs    |     Pipeline run time    |     Pipeline max memory    |
|-----------------------|--------------------------|---------------------------|-------------------|--------------------------|----------------------------|
|     Oral              |     891                  |     8.3m                  |     5417          |     2.5h                 |     34 GB                  |
|     Gut               |     192                  |     2.2m                  |     1593          |     2h                   |     30 GB                  |
|     Gut               |     192                  |     2.2m                  |     10917         |     5.5h                 |     30 GB                  |
|     Gut               |     192                  |     16.7m                 |     17293         |     13h                  |     87 GB                  |
|     Wastewater        |     33                   |     2.14m                 |     11462         |     12h                  |     47 GB                  |
|     Mock community    |     264                  |     12.8m                 |     84            |     4h                   |     44 GB                  |

## Frequently asked questions (FAQ) <a name="faq"></a>

* Can I restart the pipeline?

  Yes! The Nextflow pipeline can be resumed after interruption by adding the `-resume` option
  in the `nextflow run` command when you run the pipeline. Nextflow is intelligent enough to not rerun
  any steps if it does not need to. For example, if you want to manually provide the
  rarefaction/sampling depth after the pipeline finishes, rerun by adding
  `-resume --rarefaction_depth 5000` and only the steps that uses sampling/rarefaction depth will be rerun.
  Of course, any downstream steps will also be rerun.

- Why `cutadapt`?

  The naive Bayes classified in QIIME 2 requires the ASVs to be in the same sequence orientation.
  PacBio's CCS reads have random orientations out of the instrument, hence they need to
  to be oriented first; this can be done using either `lima` or `cutadapt`.
  Technically, `lima` has this capability too but it requires BAM input. There are many
  public datasets on SRA in FASTQ format and `lima` will not orient them. Due to the
  accuracy of HiFi reads, the performance difference between `lima` and `cutadapt` should be
  minimal in our experience.

  Without the read orientation, you will notice that the taxonomy assignments can produce
  strange results, such as archea assignment only at the highest taxonomy level.

- Many of my reads are lost in the denoise stage - what's going on?

  This can happen in extremely diverse communities such as soil where the ASVs are of very low abundance.
  In each sample, the reads supporting the ASV are very low and may not pass the DADA2 threshold to qualify
  as a cluster. In addition, DADA2 has a strict reads quality filter (maxEE parameter) that will filter
  away reads with relatively low accuracy.
  See [here](https://github.com/benjjneb/dada2/issues/841) and [here](https://github.com/benjjneb/dada2/issues/1164)
  for discussions on DADA2 algorithm and reads loss.

- I'm getting `Conda` "Safety" error indicating a corrupted package or that some
pipeline steps are not able to find specific command-line tools (e.g. qiime).

  Sometimes the conda cache can become corrupted if you run many workflows
  in parallel before the environment was created, thus causing conflicts between
  different workflow competing to create the same environment. We recommend running
  the test dataset above and then wait for it to finish first so the conda
  environment is created successfully. Subsequent runs will use the same
  environment, and will not need to recreate the environment. If the errors still
  occur, try running `conda clean -a` and remove the offending `conda` package cache
  in the cache directory. For example, if the error happens for QIIME, delete any folder in `conda info`
  "package cache" containing QIIME.

  You can try to install the QIIME 2 environment directly to inspect any error messages:

  `conda env create -n q2_test -f qiime2-2022.2-py38-linux-conda.yml`

- I've received/downloaded 16S FASTQ files that already have the primers trimmed. Can I skip primers removal?

  We recommend using the pipeline to trim the primers as it works well for HiFi sequencing
  data. However, there are many public dataset that may already have the full length
  primers trimmed, in which case you can specify `--skip_primer_trim` to skip
  primer trimming. If unsure, run with the default pipeline and the cutadapt demultiplexing rate
  (in the file `results/samples_demux_rate.tsv`) should be close to zero for all
  samples if the primers are already trimmed.

- How does the taxonomy classification work?

  The "besttax" assignment uses the `assignTaxonomy` Naive-Bayes classifier function from DADA2
  to carry out taxonomy assignment. It uses 3 databases to classify the ASVs
  (requiring a minimum bootstrap of 80 using the minBoot parameter) and the priority of assignment
  is GreenGenes2, followed by GTDB, then lastly Silva. This means, for example,
  if an ASV is not assigned at Species level using GreenGenes2, it will check if it can be assigned
  with GTDB. This ensures that we assign as many ASVs as possible.

  This process is done first at Species level, then at Genus level.
  Note that while this method will assign a high amount of ASVs, there may be issues
  such as how the taxonomy is annotated in different databases. You can change the top priority database
  by setting the `db_to_prioritize` parameter in the `nextflow.config` file. If you set it to "GTDB",
  the taxonomy assignment will prioritize GTDB -> GreenGenes2 -> Silva. If you set it to "Silva",
  the taxonomy assignment will prioritize Silva -> GreenGenes2 -> GTDB. In the `nb_tax` output folder,
  you can find the taxonomy assignment for each individual database, too, if you would like to
  use just a single database for Naive-Bayes classification.

  There is also a VSEARCH taxonomy classification using the a single database (GTDB r220 by default) only in the file called
  `results/vsearch_merged_freq_tax.tsv` that may provide a more consistent annotation. This uses
  the `classify-consensus-vsearch` plugin from `QIIME 2` and we use the "top-hits" approach
  with a stringent default hit criteria (97% identity) to classify the taxonomy of ASVs.

  The final report will contain statistics from either types of assignment. If you notice a large
  discrepancy, it can be because one method fails to assign a large amount of ASVs from the
  same genus/species. This is likely a database-related bias.

- Some species in the MSA 1003 demo data are missing!

  If you run this pipeline by default with PacBio's publicly available
  192-plex replicates ATCC-MSA1003, some 0.02% bacteria may be missing depending on
  which replicates you use due to the default `min_asv_sample` and `min_asv_totalfreq`
  parameters. These bacteria may only have a few reads in 1/2 samples, so they are
  prone to getting filtered out. You can set the two parameters to 0 to disable
  filtering and the bacteria should pop out. However, in a real dataset this
  may result in more false-positives.
  
- The percentage reads classified at species is higher than genus!
  
  You have likely bumped into strange issues with the database. For example, there are some microbes
  that have the taxonomy populated at species level, but all the other levels are empty. Unfortunately,
  database curation is out of the scope of this pipeline.

- Can I manually download the databases for taxonomic classification?

  The pipeline taxonomy classification step requires a few databases that will be downloaded with the
  `--download_db` parameters into a "databases" folder. These databases can also be downloaded
  manually from the following links if the download script above does not work. The GTDB database
  for VSEARCH will require some processing using the `QIIME 2` package. See `scripts/download_db.sh` for details.

- I want to understand more about OTU versus ASV.

  Zymo Research provides a good article on the difference between ASV and OTU here: (<https://www.zymoresearch.com/blogs/blog/microbiome-informatics-otu-vs-asv>).
  In addition, [this thread](https://forum.qiime2.org/t/esv-vs-otu-deflation-qiime1-vs-2/14867/2) on the
  `QIIME 2` forum discusses the difference in numbers through traditional OTU compared to the ASV approach.

- Can I classify with X database?
  
  With the current implementation, it is straightforward to import any database to use with VSEARCH. You will need the
  `X.fasta` sequences and the corresponding taxonomy for each sequence in the FASTA file. The taxonomy format should be
  in a 2-columns TSV file. Example:

  ```
  seq1  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Enterobacter;s__Enterobacter asburiae_B
  seq2  d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Enterobacter;s__Enterobacter asburiae_B
  ```

  Then, import both the sequences and taxonomy into the `QZA` format for `QIIME 2`:

  ```
  qiime tools import --type 'FeatureData[Sequence]' --input-path 'X.fasta' --output-path X.qza
  qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path X.taxonomy.tsv --output-path X.taxonomy.qza
  ```

  Then supply `X.qza` to `--vsearch_db`, `X.taxonomy.qza` to `--vsearch_tax`. It is not straightforward to use it with
  the Naive-Bayes approach, yet, so please also set `--skip_nb` to use only VSEARCH for classification.

- The pipeline failed in report generation with error code 137 (e.g. issue #21).

  If you have many samples with diverse sample types such as environmental samples, it is possible that DADA2 will generate
  a very large number of ASVs. The script to produce the report may subsequently fail due to running out of memory
  trying to process all the ASVs. You may want to consider splitting the different sample types into individual
  Nextflow runs to avoid this issue. Alternatively, if you have a cluster with a lot of memory, you can assign higher
  memory to the step that fails using `nextflow.config`.

- Can I get the output in hard copy instead of symlinks? (Issue #22)
  
  By default, `Nextflow` provides output in absolute symlinks (linked to files in the `work` folder) to avoid duplicating
  files. This is controlled by the `publishDir` directive (see here: <https://www.nextflow.io/docs/latest/process.html#publishdir>)
  in each process. The pipeline implements a global `--publish_dir_mode` that allows user to specify a global `publishDir`
  mode. For example, `nextflow run main.nf --publish_dir_mode copy` will provide all outputs in hard copy. Note that the files
  will still exist as duplicates in the `work` folder. You may delete the `work` folder when the pipeline finishes successfully.

- Can I change the tmpdir of QIIME 2?

  The `QIIME 2` plugin uses the `TMPDIR` environment variable to store temporary files. If you want to change the location
  of the temporary files, you can set the `TMPDIR` environment variable in the `nextflow.config` file by adding the following block:

  ```
  env {
    TMPDIR = "/path/to/tmpdir"
  }
  ```

  If you are using `singularity` or `docker`, you will have to bind the `TMPDIR` to the container. For example, if you are using
  `singularity`, you can add the following to the `nextflow.config` file:

  ```
  singularity {
    singularity.runOptions = "--bind $TMPDIR:/tmp"
  }
  ```

## References

### QIIME 2

* Bolyen, E. et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol 37, 852–857 (2019).
- For individual citations of plugins, you can use the `--citations` command for the relevant plugins. For example, if you want
  to cite `VSEARCH` plugin, type `qiime feature-classifier classify-consensus-vsearch --citations` after activating the conda
  environment. You can activate the environment installed by the pipelines by typing `conda activate $HOME/nf_conda/$ENV` (Change
  `$ENBV` to the name of the environment you want to activate).

### DADA2

* Callahan, B. J. et al. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13, 581–583 (2016).

### Seqkit

* Shen, W., Le, S., Li, Y. & Hu, F. SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11, e0163962 (2016).

### Cutadapt

* Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17, 10–12 (2011).

### GTDB database

* Parks, D. H. et al. A standardized bacterial taxonomy based on genome phylogeny substantially revises the tree of life. Nat Biotechnol 36, 996–1004 (2018).
- Parks, D. H. et al. A complete domain-to-species taxonomy for Bacteria and Archaea. Nat Biotechnol 38, 1079–1086 (2020).

### SILVA database

* Yilmaz, P. et al. The SILVA and “All-species Living Tree Project (LTP)” taxonomic frameworks. Nucleic Acids Research 42, D643–D648 (2014).
- Quast, C. et al. The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucleic Acids Research 41, D590–D596 (2013).

### RDP database

* Cole, J. R. et al. Ribosomal Database Project: data and tools for high throughput rRNA analysis. Nucleic Acids Res 42, D633–D642 (2014).

### GreenGenes2 database

* McDonald, D., Jiang, Y., Balaban, M. et al. Greengenes2 unifies microbial data in a single reference tree. Nat Biotechnol 42, 715–718 (2024). 10.1038/s41587-023-01845-1

### Krona plot

* Bd, O., Nh, B. & Am, P. Interactive metagenomic visualization in a Web browser. BMC bioinformatics 12, (2011).
- We use the QIIME 2 plugin implementation here: <https://github.com/kaanb93/q2-krona>

### Phyloseq and Tidyverse (for HTML report visualization)

* McMurdie, P. J. & Holmes, S. phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. PLOS ONE 8, e61217 (2013).
- Wickham, H. et al. Welcome to the Tidyverse. Journal of Open Source Software 4, 1686 (2019).

### PICRUSt2

* Douglas, G. M. et al. PICRUSt2 for prediction of metagenome functions. Nat Biotechnol 38, 685–688 (2020).

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
