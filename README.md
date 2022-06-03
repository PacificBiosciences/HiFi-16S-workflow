# HiFi Full-length 16S analysis with pb-16S-nf

## Workflow overview and output
This Nextflow pipeline is designed to process PacBio HiFi full-length 16S data into a set of high
quality amplicon sequence variants (ASVs) using `QIIME 2` and `DADA2`. It provides a set of visualization 
through the `QIIME 2 framework` for interactive plotting, but also generates a HTML report to report
the important statistics and top taxonomies for the 16S communities. The outputs and stages of this pipeline 
are documented [here](https://github.com/proteinosome/pb-16S-nf/blob/main/pipeline_overview.md).

We provide an example of report generated using this pipeline based on 8 replicates from the ATCC MSA-1003 
mock community sequenced on Sequel II. Right click this 
[link](https://github.com/proteinosome/pb-16S-nf/raw/main/examples/results/visualize_biom.html) and save it on 
your computer, then double click to open it. All other important outputs from the pipeline are available
in the [`examples`](https://github.com/proteinosome/pb-16S-nf/tree/main/examples) folder when you clone this repository.

## Installation and usage
This pipeline runs using Nextflow (Version 22 and above). All softwares dependencies are
managed via `Conda`. We recommend installing [`mamba`](https://github.com/mamba-org/mamba)
to speed up the conda environment installation. The default `nextflow.config` file 
enables the use of `mamba` by default. You can install Nextflow following the instruction
from Nextflow [documentation](https://www.nextflow.io/docs/latest/getstarted.html) or via Conda:

```
# (Optional but recommended) Install mamba
conda install mamba -n base -c conda-forge
conda install -c bioconda nextflow
```

After installing Nextflow and `mamba`(Optional but recommended), clone the repository and
download databases using the following commands:

```
git clone https://github.com/proteinosome/pb-16S-nf.git
cd pb-16S-nf
bash scripts/download_db.sh
```

The taxonomy classification step of the pipeline requires a few databases that will be downloaded with the
`download_db.sh` script above into a `references` folder. These databases can also be downloaded
manually from the following links if the download script above does not work.

- `--vsearch_db` and `--vsearch_tax` provided by the `QIIME 2` community
  - [`silva-138-99-seqs.qza`](https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza)
  - [`silva-138-99-tax.qza`](https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza)
- `--silva_db` provided by `DADA2`
  - [`silva_nr99_v138.1_wSpecies_train_set.fa.gz`](https://zenodo.org/record/4587955)
- `--gtdb_db` provided by `DADA2`
  - [`GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz`](https://zenodo.org/record/4735821)
- `--refseq_db` provided by `DADA2`
  - [`RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz`](https://zenodo.org/record/4735821)

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
  
  nextflow run main.nf --input samples.tsv --metadata metadata.tsv \\
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
                 (default: 100)
  --min_asv_totalfreq    Total frequency of any ASV must be above this threshold
                         across all samples to be retained. Set this to 0 to disable filtering
                         (default 5)
  --min_asv_sample    ASV must exist in at least min_asv_sample to be retained. 
                      Set this to 0 to disable. (default 2)
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
  --vsearch_db	Location of VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Location of VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --silva_db   Location of Silva 138 database for taxonomy classification 
  --gtdb_db    Location of GTDB r202 for taxonomy classification
  --refseq_db    Location of RefSeq+RDP database for taxonomy classification
  --skip_primer_trim    Skip all primers trimming (switch off cutadapt and DADA2 primers
                        removal) (default: trim with cutadapt)
  --colorby    Columns in metadata TSV file to use for coloring the MDS plot
               in HTML report (default: condition)
```

To test the pipeline, run this example below. Note that the path of the database needs
to be changed to their respective locations on your server if it's different (See parameters above). If you 
follow the command above, the databases will be downloaded into a `references` folder in the `pb-16S-nf` folder
and you do not need to specify the path.

Please change the path on "test_sample.tsv" to point to the correct absolute location of the test dataset
(You can find it by typing `readlink -f test_data/test_1000_reads.fastq.gz` when you're in
the cloned pb-16S-nf folder). Conda environment will by default be created at 
`$HOME/nf_conda` folder unless changed in the `nextflow.config` file. Once the conda environment
is created it will be reused by any future run.

```
# Find the path of test dataset and replace the path in test_sample.tsv
readlink -f test_data/test_1000_reads.fastq.gz

nextflow run main.nf --input test_sample.tsv \
    --metadata test_metadata.tsv -profile conda \
    --dada2_cpu 32 --vsearch_cpu 32 --outdir results \
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
required. See the `report_results` folder created after the pipeline finishes running 
for DAG and resources report.

## Frequently asked questions (FAQ)
* Can I restart the pipeline?

  Yes! Nextflow pipeline can be resumed after interruption by adding the `-resume` option 
  in the `nextflow run` command when you run the pipeline. Nextflow is smart enough to not rerun
  any step if it does not need to. For example, if you want to manually provide the 
  rarefaction/sampling depth after the pipeline finishes, rerun by adding 
  `-resume --rarefaction_depth 5000` and only the steps that uses sampling/rarefaction depth will rerun.
  Of course, any step downstream of those will also rerun.

* Why `cutadapt`?

  The naive Bayes classified in QIIME 2 requires the ASVs to be in the same sequence orientation.
  PacBio's CCS reads have random orientations out of the instrument, hence they need to
  to be oriented first, and this can be done with either `lima` or `cutadapt`. 
  Technically, `lima` has this capability too but it requires BAM input. There are many
  public dataset on SRA that is in FASTQ format and `lima` will not orientate them. Due to the 
  accuracy of HiFi reads, the performance difference between `lima` and `cutadapt` should be
  minimal in our experience.

  Without the read orientation, you will notice that the taxonomy assignments can produce
  weird results (e.g. archea assignment only at the highest taxonomy level).

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
  This pipeline uses 3 databases to classify the ASVs (requiring a minimum bootstrap of 80 using
  the minBoot parameter) and the priority of assignment
  is Silva 138, followed by GTDB r202, then lastly RefSeq + RDP. This means for example
  if an ASV is not assigned at Species level using Silva, it will check if it can be assigned
  with GTDB. This ensure we assign as many ASVs as possible. 

  This process is done first at Species level, then at Genus level. In addition, if any ASV
  is assigned as "uncultured" or "metagenome", it will go through the iterative assignment
  process just like the unclassified ASVs.

  There is also a VSEARCH taxonomy classification using Silva database only in the file called 
  `results/vsearch_merged_freq_tax.tsv` that may work better in some cases. The final 
  report will contain statistics from either types of assignment. If you notice a large
  discrepancy, it can be because one method fail to assign a large amount of ASVs from the
  same genus/species. This is likely a database-related bias.

* Some species in MSA 1003 demo data are missing!

  If you run this pipeline by default with PacBio's publicly available
  192-plex replicates ATCC-MSA1003, some 0.02% bacteria may be missing depending on
  which replicates you use due to the default `min_asv_sample` and `min_asv_totalfreq`
  parameters. These bacteria may only have a few reads in 1/2 samples, so they're
  prone to getting filtered out. You can set the two parameters to 0 to disable
  filtering and the bacteria should pop out. However, in real dataset this 
  may result in more false-positives.

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
