# pb-16S-analysis Pipeline using QIIME 2

This pipeline requires java 11 (Can be installed via Conda) to use Nextflow.
Nextflow version 22 onwards is needed as the pipeline is written in DSL 2 language.

QIIME 2 installation via conda using yml file. 

* How to run with trace reports and visualization of workflow etc:
`nextflow run ~/git/pacbio/pb-16S-nf/main.nf --input sample.tsv \
    --metadata metadata.tsv -profile conda \
    --rarefaction_depth 10000 \
    -with-trace -with-timeline timeline.html -with-dag workflow.html \
    --dada2_cpu 32 --vsearch_cpu 32 \
    -with-report run_report.html`
    
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
