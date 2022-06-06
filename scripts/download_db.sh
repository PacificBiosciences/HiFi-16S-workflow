#!/bin/bash

echo "Downloading SILVA sequences for VSEARCH..."
wget -N -P references/ --content-disposition 'https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1'
echo "Downloading GTDB sequences for VSEARCH..."
wget -N -P references/ --content-disposition 'https://zenodo.org/record/4735821/files/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz?download=1'
echo "Downloading RefSeq + RDP sequences for VSEARCH..."
wget -N -P references/ --content-disposition 'https://zenodo.org/record/4735821/files/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz?download=1'

echo "Downloading SILVA sequences and taxonomies for Naive Bayes classifier"
wget -N -P references/ --content-disposition 'https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza'
wget -N -P references/ --content-disposition 'https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza'

echo "Downloading GTDB database"
# wget -N -P references/ --content-disposition 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_all/ssu_all_r207.tar.gz'
# tar xfz references/ssu_all_r207.tar.gz
# mv ssu_all_r207.fna references/ssu_all_r207.fna
# grep "^>" references/ssu_all_r207.fna | awk '{print gensub(/^>(.*)/, "\\1", "g", $1),gensub(/^>.*\ (d__.*)\ \[.*\[.*\[.*/, "\\1", "g", $0)}' OFS=$'\t' > references/ssu_all_r207.taxonomy.tsv
wget -N -P references/ --content-disposition 'https://zenodo.org/record/6616062/files/GTDB_ssu_all_r207.qza?download=1'
wget -N -P references/ --content-disposition 'https://zenodo.org/record/6616062/files/GTDB_ssu_all_r207.taxonomy.qza?download=1'
