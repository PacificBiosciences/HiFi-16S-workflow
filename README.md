# pb-16S-analysis
Information and tools for the analysis of fully length 16S data. 

* Updated 2022-3-26: First upload
* Updated 2022-4-13: Example instruction to use other databases

- [Analyzing PacBio HiFi Mock Community 16S Data with `QIIME 2`](#analyzing-pacbio-hifi-mock-community-16s-data-with-qiime-2)
  - [Introduction](#introduction)
  - [ATCC MSA-1003 Mock Community](#atcc-msa-1003-mock-community)
  - [Importing Data into QIIME 2](#importing-data-into-qiime-2)
  - [Using `QIIME 2` View](#using-qiime-2-view)
  - [Denoising HiFi Reads into ASVs](#denoising-hifi-reads-into-asvs)
  - [Summarizing Denoised Statistics](#summarizing-denoised-statistics)
  - [Rarefaction plot](#rarefaction-plot)
  - [Classifying Taxonomy of ASVs](#classifying-taxonomy-of-asvs)
  - [Taxonomy bar plot](#taxonomy-bar-plot)
  - [Phylogeny analysis](#phylogeny-analysis)
  - [Diversity analysis and PCoA](#diversity-analysis-and-pcoa)

## Introduction
`QIIME 2` is traditionally the tool of choice for 16S researcher due to the large amount of plugins and its ease of use. Almost everything that one will need to analyze complex communites based on 16S rRNA sequencing can be done within the QIIME 2 platform. However, with PacBio full-length 16S, the analysis workflow typically requires one to use tool such as `DADA2` in R first, then export the results into QIIME 2 format.

Thankfully, the recent release of `QIIME 2` version 2022-2 has now integrated a wrapper for DADA2 that allows one to process CCS data directly without the need of using R (Although behind the scene, it's just running DADA2 which is based on R!). For this tutorial, we will use the public demo dataset ATCC-MSA1003 16S sequencing provided on PacBio's website to demonstrate a typical analysis workflow taking FASTQ files to taxonomy barplot.

`QIIME 2` [documentation page](<https://docs.qiime2.org/>) provides many useful and comprehensive tutorials that we highly encourage any 16S users to go through. The tutorial we provide here is created as an example of how you can reconstruct the mock community dataset using a publicly available HiFi 16S dataset.

***

## ATCC MSA-1003 Mock Community
A mock community ATCC-MSA1003 has been sequenced on a Sequel II system. The same sample is replicated 192 times to simulate a 192-plex multiplexed experiment, commonly done for 16S rRNA sequencing due to the high throughput of Sequel II. This mock community contains 20 strains of bacteria mixed in a staggered fashion, where the bacteria abundance ranges from 0.02% to 18%. For details on the mock community, please refer to ATCC website: <https://www.atcc.org/products/msa-1003>.

The full 192-plex dataset can be downloaded here: [DevNet](<https://github.com/PacificBiosciences/DevNet/wiki/16S-Data-Set-Sequel-II-System-2.0-Release>). We will use 8 of the 192 barcodes pairs to save on the computational time for this demo. We will use these 8 FASTQ files for the tutorial and results below:
```
demultiplex.16S_For_bc1005--16S_Rev_bc1056.hifi_reads.fastq.gz
demultiplex.16S_For_bc1005--16S_Rev_bc1057.hifi_reads.fastq.gz
demultiplex.16S_For_bc1005--16S_Rev_bc1062.hifi_reads.fastq.gz
demultiplex.16S_For_bc1005--16S_Rev_bc1075.hifi_reads.fastq.gz
demultiplex.16S_For_bc1005--16S_Rev_bc1100.hifi_reads.fastq.gz
demultiplex.16S_For_bc1007--16S_Rev_bc1075.hifi_reads.fastq.gz
demultiplex.16S_For_bc1020--16S_Rev_bc1059.hifi_reads.fastq.gz
demultiplex.16S_For_bc1024--16S_Rev_bc1111.hifi_reads.fastq.gz
```
where each FASTQ file was the results of demultiplexing the full 192-plex dataset. Create a 8-plex folder and put the FASTQ files inside.

* Additional homework: What is the average quality of the sequences? Can you use tools such as `seqkit` and `datamash` to find out?

***

## Importing Data into QIIME 2
`QIIME 2` is straightforward to install it via Anaconda. Please refer to `QIIME 2` webpage [here](<https://docs.qiime2.org/2022.2/install/>) for installation instructions. After installing `QIIME 2`, activate the environment.

To test if `QIIME 2` has been installed, type `qiime --help` and you will see a block of texts documenting the usage of `QIIME 2` and the various functions or plugins available. The first tool we will use is `qiime tools import` that allows us to import the FASTQ file into a `QIIME 2 qza` format. 

In order to import the data, we will need a sample "manifest" file that allows us to give each of the FASTQ file a sample name. For example, to name the file according to the barcodes pair in the file names:
```
echo -e 'sample-id\tabsolute-filepath' > sample.manifest

for i in 8-plex/*.fastq.gz; do fname=$(basename ${i}); sample=$(echo ${fname} |\
 awk '{print gensub(/.*_(bc.*)--.*(bc.*).hifi_reads.fastq.gz/, "\\1_\\2", "g", $0)}');\
 echo -e "${sample}\t$(readlink -f ${i})" >> sample.manifest;\
 done
```

You can open the sample manifest file to look at how it is formatted. It's a simple two-columns TSV file giving a sample name for each of the FASTQ file that we are going to analyze. To import the FASTQ file, we can type:

```
qiime tools import --type 'SampleData[SequencesWithQuality]' \
 --input-path sample.manifest \
 --output-path samples.qza \
 --input-format SingleEndFastqManifestPhred33V2
```
which tells `QIIME 2` to import FASTQ file using the sample manifest file into a file called `samples.qza`. It also specifies the import file format and type. This should take only a few seconds to complete. Once completed, you will see an output saying:
```
Imported sample.manifest as SingleEndFastqManifestPhred33V2 to samples.qza
```
and if you type `ls` now you will see a new file created called `samples.qza`. This file can now be used for all the downstream analysis. You can summarize the read statistics of each barcode pair by running:
```
qiime demux summarize --i-data samples.qza \
 --o-visualization samples.demux.summary.qzv
```
which again, should finish very quickly (< 1min). This will create a file called `samples.demux.summary.qzv` that you can download to your local computer, e.g.(on your own computer):
```
scp user@17.xx.xx.xx:/path/to/samples.demux.summary.qzv .
```

Usually, in an 16S experiment, one will also have multiple conditions for different groups of samples. `QIIME 2` can use a file specifying the samples conditions for downstream analysis and visualization later on. Let's create two artificial groups here where one group is called RepA and one is called RepB. The samples are all identical, so this is just artificial grouping:

```
awk 'FNR==1{print "sample_name","condition"}; FNR>1 && FNR <6 {print $1,"RepA"}; FNR>5 {print $1, "RepB"}' OFS=$'\t' sample.manifest > metadata.tsv
```

***

## Using `QIIME 2` View
One of the nice things about `QIIME 2` is that it comes with a web browser that allows you to visualize many of the plots/statistics produced by the `QIIME 2` package (typically ending with a `.qzv` extension, similar to what we have just generated above for the reads statistics). You can use your web browser and navigate to <https://view.qiime2.org/>. Next, drag and drop the file `samples.demux.summary.qzv` into the browser and you can now see the read statistics for the dataset we've imported (note that with PacBio data there's no forward and reverse sequence, so all reads will be seen under "forward sequence"). Here, each barcode pair has on average 13k HiFi reads. 

***

## Denoising HiFi Reads into ASVs
As discussed above, the current recommendations for analyzing full-length 16S HiFi reads is to use `DADA2` (Original [paper](<https://www.nature.com/articles/nmeth.3869>) which is a tool to infer high resolution amplicon sequence variants (ASV) from the reads. It is able to resolve even down to 1-2 nucleotides difference, and hence is able to take advantage of the high accuracy of HiFi reads. For example, you can run `DADA2` using R script following the tutorial provided by `DADA2` author, e.g. On the [Zymo mock community](<https://benjjneb.github.io/LRASManuscript/LRASms_Zymo.html>).

The latest version of `QIIME 2` has now incorporated a function called `dada2 denoise-ccs` that can be used to process HiFi reads. Let's use this command to infer ASVs from our HiFi reads (This should take ~ 5 to 8 mins with 8 threads):

```
qiime dada2 denoise-ccs --i-demultiplexed-seqs ./samples.qza \
 --o-table dada2-ccs_table.qza \
 --o-representative-sequences dada2-ccs_rep.qza \
 --o-denoising-stats dada2-ccs_stats.qza \
 --p-min-len 1000 --p-max-len 1600 \
 --p-front AGRGTTYGATYMTGGCTCAG --p-adapter RGYTACCTTGTTACGACTT \
 --p-n-threads 8
```

A few notable option for `dada2 denoise-ccs` (Type `qiime dada2 denoise-ccs` for an exhaustive list):

1. `--p-min-len 1000` and `--p-max-len 1600` specifies the minimum and maximum sequence length. For 16S sequences, `QIIME 2` suggests using 1000 and 1600, respectively. 
2. `--p-front AGRGTTYGATYMTGGCTCAG` and `--p-adapter RGYTACCTTGTTACGACTT` which are the full length 16S forward (F27) and reverse (R1492) primers. This will allow DADA2 to trim away the primer region.
3. `--p-n-threads 8` specifies that the software can use up to 8 threads.

***

## Summarizing Denoised Statistics
After denoising, it is important to check the read statistics to see how many reads have been filtered and what is the amount of chimeric reads (even after filtering). `QIIME 2` can do this using:

```
qiime metadata tabulate --m-input-file ./dada2-ccs_stats.qza \
 --o-visualization ./dada2-ccs_stats.qzv
```

Again, you can use `QIIME 2 View` to view the table (filename `dada2-ccs_stats.qzv`) in your web browser. You should observe around 70% of reads categorized as non-chimeric reads. This number can vary from sample to sample.

| sample-id      | input | primer-removed | percentage of input primer-removed | filtered | percentage of input passed filter | denoised | non-chimeric | percentage of input non-chimeric |
| -------------- | ----- | -------------- | ---------------------------------- | -------- | --------------------------------- | -------- | ------------ | -------------------------------- |
| bc1005\_bc1056 | 13682 | 12017          | 87.83                              | 9823     | 71.8                              | 9746     | 9746         | 71.23                            |
| bc1005\_bc1057 | 14009 | 12363          | 88.25                              | 10128    | 72.3                              | 10042    | 10040        | 71.67                            |
| bc1005\_bc1062 | 13051 | 11523          | 88.29                              | 9561     | 73.26                             | 9499     | 9497         | 72.77                            |
| bc1005\_bc1075 | 13626 | 12003          | 88.09                              | 9825     | 72.1                              | 9740     | 9740         | 71.48                            |
| bc1005\_bc1100 | 12494 | 11019          | 88.19                              | 9097     | 72.81                             | 9025     | 9021         | 72.2                             |
| bc1007\_bc1075 | 14072 | 12497          | 88.81                              | 10190    | 72.41                             | 10098    | 10098        | 71.76                            |
| bc1020\_bc1059 | 13717 | 12050          | 87.85                              | 9860     | 71.88                             | 9801     | 9801         | 71.45                            |
| bc1024\_bc1111 | 12887 | 11315          | 87.8                               | 9301     | 72.17                             | 9215     | 9213         | 71.49                            |

Another tables that summarize the features (ASVs here) observed for each sample can be generated using:

```
qiime feature-table summarize --i-table ./dada2-ccs_table.qza \
 --m-sample-metadata-file metadata.tsv \
 --o-visualization ./dada2-ccs_table.qzv
```

which also makes use of the metadata file we generated above. A useful table here would be the table that summarizes the number of reads assigned to each feature. 

|                   | Frequency |
| ----------------- | --------- |
| Minimum frequency | 4         |
| 1st quartile      | 84        |
| Median frequency  | 147       |
| 3rd quartile      | 2,201.00  |
| Maximum frequency | 11,733.00 |
| Mean frequency    | 1,455.77  |

The second tab of this visualization allows you to interactively select a cut-off where you can see how many samples "drop-out" at a certain cut-off. Finally, you can find a table on the third tab that provides detailed information on how many samples are the ASVs detected in, and how many reads in total for each of these ASVs. For example, the first 5 ASVs here are:

|                                  | Frequency | # of Samples Observed In |
| -------------------------------- | --------- | ------------------------- |
| b80a1bd6b19379171cafb03a4ddb26a3 | 11,733    | 8                         |
| 61de72098cbf9f693f09c863d1c5a586 | 11,111    | 8                         |
| a0530711bdec8d611872863731f96dde | 8,723     | 8                         |
| bc7d5e82dae59e2e5e8321bb5c5eae4b | 5,208     | 8                         |
| 13bc41e35faa4f40a01de5a972c0c4ac | 4,936     | 8                         |
| f3cfb383dcfb753cad226072f6665784 | 4,400     | 8                         |
| 7dcace79bc8188077e07b0085aa57b7f | 3,824     | 8                         |
| 09a8fa70427aa68f16934c729b826071 | 2,459     | 8                         |

***

## Rarefaction plot
One of the very first analysis for 16S experiment is to plot a rarefaction curve. This would allow us to understand if we're saturating the amount of features (or ASVs) that can be discovered from the sample. Let's do that with a simple command here:

```
qiime diversity alpha-rarefaction --i-table dada2-ccs_table.qza \
 --m-metadata-file metadata.tsv \
 --o-visualization ./alpha-rarefaction-curves.qzv \
 --p-min-depth 10 --p-max-depth 10000
```

<img src="https://github.com/PacificBiosciences/pb-16S-analysis/blob/main/figures_16S_tutorial/rarefaction.png" alt="alpha_rarefaction" width="1440"/>

At what sequencing depth do you think we are saturating the experiment? It is important to note that we have a very simple mock community here. In real samples, the diversity and number of ASVs are often much higher than what we observe here. In addition, there's another plot that allows us to check for the number of samples that drop-out beyond a certain depth (useful for QC if you have a large amount of samples!).

***

## Classifying Taxonomy of ASVs
A popular choice of classifier for 16S sequences is the naive Bayesian based `qiime feature-classifier classify-sklearn`. However, in our experience, this is not very sensitive for our full-length mock community ASVs. `DADA2` comes with a tool called `assignTaxonomy` that uses a similar naive Bayesian classifier approach which works well (in the mock community), but requires one to use `DADA2` in R. 

With `QIIME 2`, we are going to use `VSEARCH` to classify the taxonomy of each of the ASVs from DADA2. `VSEARCH` is a significant faster alternative to `BLAST` and is able to classify almost all of the species in our mock community here. We will use the `SILVA` version 138 database here for classification. It is debatable which database is the best, as they can all perform differently depending on the communities. We chose `SILVA` here as it's continuously updated. The database formatted for `QIIME 2` can be downloaded from [`QIIME 2` resources](<https://docs.qiime2.org/2022.2/data-resources/>) webpage.

```
qiime feature-classifier classify-consensus-vsearch --o-classification ./taxonomy.vsearch.qza \
 --i-query dada2-ccs_rep.qza \
 --i-reference-reads references/silva-138-99-seqs.qza \
 --i-reference-taxonomy references/silva-138-99-tax.qza  \
 --p-threads 8 \
 --p-maxrejects 100 --p-maxaccepts 10 --p-perc-identity 0.97
```

This should take < 1 min with 8 threads. Three important parameters that we changed for the classifier are `--p-maxrejects 100 --p-maxaccepts 10 --p-perc-identity 0.97`. `--p-max-rejects 100` specify that the search for taxonomy hit will stop after 100 non-matching sequences, and `--p-maxaccepts 10` will stop the search if there's 10 taxonomy hits. In combination with `--p-min-consensus` which has a default value of 0.51, more than half of the "hits" (set to be 10 here) need to have the same assignment at any specific taxonomy level for the sequences to be classified confidently. Finally, `--p-perc-identity 0.97` specifies that the taxonomy hits need to be at least 97% identical. These options ensure that `VSEARCH` runs very fast (Certainly much faster than the traditional `BLAST` approach) and keeps only highly identical hits.

Let's look at what are the taxonomies assigned. You can use the first command to export it into the usual `qzv` form for `QIIME 2 View`, or you can use the second command to export it into a TSV format (useful for programmatic manipulation with R/Python etc):

```
qiime metadata tabulate --m-input-file ./taxonomy.vsearch.qza \
 --o-visualization taxonomy.vsearch.qzv

qiime tools export --input-path taxonomy.vsearch.qza \
 --output-path tax_export
```

The top 5 rows in `QIIME 2 View` for example showed the ASVs IDs and the taxonomy of the ASV, in addition to the confidence of the assignment based on the number of top 10 hits (discussed above). As you can see, some ASVs cannot be assigned to species level and will end at the genus level for its taxonomic classification.

| Feature ID                       | Taxon                                                                                                                                                           | Consensus   |
| -------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------- |
| #q2:types                        | categorical                                                                                                                                                     | categorical |
| 021408f68ac3d943f7e16261373b2a61 | d\_\_Bacteria; p\_\_Firmicutes; c\_\_Clostridia; o\_\_Clostridiales; f\_\_Clostridiaceae; g\_\_Clostridium\_sensu\_stricto\_1                                   | 1           |
| 04d3de55f475a0888205c9507901a096 | d\_\_Bacteria; p\_\_Actinobacteriota; c\_\_Actinobacteria; o\_\_Actinomycetales; f\_\_Actinomycetaceae; g\_\_Actinomyces                                        | 1           |
| 070e52443d30c5d9e99c8cb06412db1a | d\_\_Bacteria; p\_\_Firmicutes; c\_\_Clostridia; o\_\_Clostridiales; f\_\_Clostridiaceae; g\_\_Clostridium\_sensu\_stricto\_1; s\_\_Clostridium\_beijerinckii   | 0.6         |
| 07b2978940b50886f3866b04b7f982cf | d\_\_Bacteria; p\_\_Firmicutes; c\_\_Bacilli; o\_\_Lactobacillales; f\_\_Streptococcaceae; g\_\_Streptococcus; s\_\_Streptococcus\_agalactiae                   | 1           |
| 09a8fa70427aa68f16934c729b826071 | d\_\_Bacteria; p\_\_Proteobacteria; c\_\_Gammaproteobacteria; o\_\_Enterobacterales; f\_\_Enterobacteriaceae; g\_\_Escherichia-Shigella; s\_\_Escherichia\_coli | 0.7         |
| 0bf3e0641c88ec588bc66b567c3fa63d | d\_\_Bacteria; p\_\_Proteobacteria; c\_\_Gammaproteobacteria; o\_\_Enterobacterales; f\_\_Enterobacteriaceae; g\_\_Escherichia-Shigella; s\_\_Escherichia\_coli | 0.8         |
| 0ccaf57fc10e3ca943b9911e91c05de1 | d\_\_Bacteria; p\_\_Firmicutes; c\_\_Clostridia; o\_\_Clostridiales; f\_\_Clostridiaceae; g\_\_Clostridium\_sensu\_stricto\_1; s\_\_Clostridium\_beijerinckii   | 0.8         |
| 0ed2d53fd823b2ef922b6f349cade1cd | d\_\_Bacteria; p\_\_Proteobacteria; c\_\_Gammaproteobacteria; o\_\_Burkholderiales; f\_\_Neisseriaceae; g\_\_Neisseria; s\_\_Neisseria\_meningitidis            | 1           |

* Additional knowledge: If you use DADA2 to classify, the default recommendation for species assignment is to use `assignSpecies` which requires perfect match between the ASVs and the reference sequences. While this is important for shorter amplicons such as V3-V4 to reduce false-positives, it may not work well with full-length 16S sequencing. See this GitHub issue for example: <https://github.com/benjjneb/dada2/issues/990>
* Tips: `qiime feature-table filter-features` can be used to filter out features (ASVs) that have extremely low abundance (e.g. singletons)
* Bonus subject: If you wish to use other database such as GTDB for classification, this is an example of importing the GTDB genome FASTA and taxonomy file into `qza` format to use with `VSEARCH`. Note that the instruction below was tested with GTDB r202, so it may require changes if you wish to use other versions of the database.

``` bash
# Here we will use the all genomes FASTA file from GTDB. You can also use
# the GTDB representative sequences but the amount of sequences are already
# collapsed in that set of data, so the consensus method of VSEARCH
# classification may have low "confidence" due to low number of hits and
# will require you to optimize the parameter.

# Extract the taxonomy sequences from the GTDB FASTA header using awk. Take
# a look at how the taxonomies are formatted and you can similarly
# generate such a taxonomy file for any database.
grep "^>" ssu_all_r202.fna | \
    awk '{print gensub(/^>(.*)/, "\\1", "g", $1),gensub(/^>.*\ (d__.*)\ \[.*\[.*\[.*/, "\\1", "g", $0)}' OFS=$'\t' > ssu_all_r202.taxonomy.tsv

# Import the taxonomy
qiime tools import  --type 'FeatureData[Taxonomy]' \
    --input-path ssu_all_r202.taxonomy.tsv \
    --output-path ssu_all_r202.taxonomy.qza \
    --input-format HeaderlessTSVTaxonomyFormat

# Import the genome FASTA
qiime tools import --type 'FeatureData[Sequence]' \
    --input-path ssu_all_r202.fna \
    --output-path ssu_all_r202.qza

# Classify with VSEARCH. You may have to play around with the parameters to find out what
# works best for your sample
qiime feature-classifier classify-consensus-vsearch \
    --o-classification ./taxonomy.vsearch.qza  \
    --i-query ../dada2-ccs_rep.qza \
    --i-reference-reads bac120_ssu_reps_r202.qza \
    --i-reference-taxonomy bac120_taxonomy.qza \
    --p-threads 8 \
    --p-maxrejects 100 --p-maxaccepts 20 \
    --p-perc-identity 0.97
```

***

## Taxonomy bar plot
Another one of the most common plots we can see in 16S studies is the taxonomy bar plot. We can generate this using the `QIIME 2` command:

```
qiime taxa barplot --i-table dada2-ccs_table.qza \
 --i-taxonomy taxonomy.vsearch.qza \
 --m-metadata-file ./metadata.tsv \
 --o-visualization ./taxa_barplot.qzv
```

Visualize this in `QIIME 2 View`, and tinker around with the functions available to look at the different taxonomies in your sample. Did we find all the taxonomies in ATCC-MSA1003? How many ASVs were assigned to species level?

<img src="https://github.com/PacificBiosciences/pb-16S-analysis/blob/main/figures_16S_tutorial/qiime2_taxa_barplot.png" alt="QIIME2_taxa_bplot" width="1440"/>

To evaluate how well the results compared to the ground truth in a more systematic way, you can make use of the results from `QIIME 2` using your favourite data analysis language. Here's a plot (Plotted using `R` with the results here) showing the abundance comparing to the true abundance from ATCC-MSA1003. Note that the abundance for all 8 samples are averaged for this plot. From the `QIIME 2 View` barplots visualization, you should also see that from sample to sample, the variation in abundance is minimal.

<img src="https://github.com/PacificBiosciences/pb-16S-analysis/blob/main/figures_16S_tutorial/species_compare_truth_barplot_corr.png" alt="Compare to real abundance" width="1440"/>

It looks like we did really well except for two missing bacteria, *Bifidobacterium adolescentis* and S*chaalia odolyntica*. At 0.02%, there's only on average 2 reads out of 10k reads that come from those species. Coupled with uneven sequencing of different species (but not uncommon due to the complexity of microbes in a complex community), it is not surprising to miss them at this sequencing depth. In fact, if you combine the full 192-plex datasets, you will be able to detect these two species. The plot also showed extremely good correlation of abundance calculated using the number of CCS reads assigned by `DADA2` to each ASV compared to the true abundance across different range of abundances.

***

## Phylogeny analysis
The power of `QIIME 2`, as mentioned, is that it can carry out many useful analysis routinely done in 16S experiment. For example, we can construct a phylogenetic tree of all the ASVs generated by using:

```
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences dada2-ccs_rep.qza \
 --output-dir phylogeny-align-to-tree-mafft-fasttree
```

Note that by default, `qiime phylogeny align-to-tree-mafft-fasttree` will create both an unrooted and a rooted (midpoint rooting) tree. The `.qza` tree format can be visualized directly in iToL website: <https://itol.embl.de/upload.cgi>. After uploading the phylogenetic tree in `.qza` format, you can also upload the taxonomy file from `QIIME 2` generated from above (`taxonomy.vsearch.qza`) to annotate the tree with the taxonomy (top right corner, `Datasets` followed by `Upload annotation`). Here's a tree generated with our data:

<img src="https://github.com/PacificBiosciences/pb-16S-analysis/blob/main/figures_16S_tutorial/iToL.png" alt="iToL" width="1440"/>

***

## Diversity analysis and PCoA
Often, we want to understand the diversity of the community we're sequencing and alpha/beta diversity are two common approaches. In addition, a PCoA (Principal coordinates analysis) based on these diversity metrics can be used for visualization purpose. Lucky for us, `QIIME 2` provides a convenient one-liner to generate all we need:

```
qiime diversity core-metrics-phylogenetic --i-table ./dada2-ccs_table.qza \
 --i-phylogeny ./phylogeny-align-to-tree-mafft-fasttree/rooted_tree.qza \
 --m-metadata-file ./metadata.tsv \
 --p-sampling-depth 9000 \
 --output-dir ./core-metrics-results
```

This command will create a folder called `core-metrics-results` containing various metrics of diversity. The ["Moving Pictures" tutorial](<https://docs.qiime2.org/2022.2/tutorials/moving-pictures/>) from `QIIME 2` is a great resource for understanding what each of these metrics represent. We know our samples are all replicates of each other and we do not expect the groups that we created in the "metadata.tsv" file to be meaningful. Let's test this by comparing the beta-diversity between the two groups using Bray-Curtis distance:

```
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
 --m-metadata-file metadata.tsv --m-metadata-column condition \
 --o-visualization core-metrics-results/bray_curtis_bet_div_comp.qzv \
 --p-pairwise
```

The `--p-pairwise` parameter asks `QIIME 2` to compare between all groups (in our case there's only two groups) in the column `condition` in our metadata TSV file. The results show a p-value that's around 0.08. Purely by chance, it appears that our random label almost resulted into two groups that's a little different! (If you look at the rarefaction curve just now, one group appears to resolve into more ASVs). 

|                        | PERMANOVA results |
| ---------------------- | ----------------- |
| method name            | PERMANOVA         |
| test statistic name    | pseudo-F          |
| sample size            | 8                 |
| number of groups       | 2                 |
| test statistic         | 1.60154           |
| p-value                | 0.079             |
| number of permutations | 999               |

Question: How did we determine `--p-sampling-depth`?


## DISCLAIMER
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
