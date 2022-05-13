# pb-16S-analysis Pipeline using QIIME 2

## Installation and usage
This pipeline requires java 11 (Can be installed via Conda) to use Nextflow.
Nextflow version 22 onwards is needed as the pipeline is written in DSL 2 language.

QIIME 2 installation via conda using yml file. 

This repository contains a set of small test data as well as example samples and metadata
TSV file that can be used to test the pipeline. Note that you need to change the
path of the FASTQ in `test_samples.tsv` to point to the correct location of the
test dataset.

After installing Nextflow 22, run the following command to see the options for
the pipeline:

```
nextflow run main.nf --help
```

The taxonomy classification step of the pipeline requires a database. We recommend
using the Silva 138 database that can be downloaded at:

[silva-138-99-seqs.qza](https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza)
[silva-138-99-tax.qza](https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza)

To test the pipeline, run this example below. Note that the path of the database needs
to be changed to the location on your server.

```
nextflow run main.nf --input test_sample.tsv \
    --metadata test_metadata.tsv -profile conda \
    --rarefaction_depth 50 \
    --dada2_cpu 32 --vsearch_cpu 32 --outdir results \
    --vsearch_db /path/to/silva-138-99-seqs.qza \
    --vsearch_tax /path/to/silva-138-99-tax.qza
```


Pipeline is still under development. The nextflow.config file by default will generate workflow DAG and resources
report, so there's no need to specify on command line.:

## Frequently asked questions (FAQ)
* A lot of my reads are lost in the denoise stage, what's going on?
This can happen in extremely diverse community such as soil where the ASVs are of very low abundance.
In each sample, the reads supporting the ASV are very low and may not pass DADA2 threshold to qualify
as a cluster. See [here](https://github.com/benjjneb/dada2/issues/841) for a discussion on
the algorithm.

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
