# pb-16S-analysis Pipeline using QIIME 2

This pipeline requires java 11 (Can be installed via Conda) to use Nextflow.
Nextflow version 22 onwards is needed as the pipeline is written in DSL 2 language.

QIIME 2 installation via conda using yml file. 

This repository contains a set of small test data as well as example samples and metadata
TSV file that can be used to test the pipeline. Note that you need to change the
path of the FASTQ in `test_samples.tsv` to point to the correct location of the
test dataset.

An example of how to run with trace reports and visualization of workflow etc (
the nextflow.config file by default will generate workflow DAG and resources
report, so there's no need to specify on command line):

```
nextflow run ~/git/pacbio/pb-16S-nf/main.nf --input sample.tsv \
    --metadata metadata.tsv -profile conda \
    --rarefaction_depth 10000 \
    --dada2_cpu 32 --vsearch_cpu 32
```


Pipeline is still under development.

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
