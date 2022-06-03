#!/bin/bash

echo "Downloading SILVA sequences for VSEARCH..."
wget --content-disposition 'https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1'
echo "Downloading GTDB sequences for VSEARCH..."
wget --content-disposition 'https://zenodo.org/record/4735821/files/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz?download=1'
echo "Downloading RefSeq + RDP sequences for VSEARCH..."
wget --content-disposition 'https://zenodo.org/record/4735821/files/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz?download=1'

echo "Downloading SILVA sequences and taxonomies for Naive Bayes classifier"
wget --content-disposition 'https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza'
wget --content-disposition 'https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza'
