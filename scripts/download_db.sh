#!/bin/bash

echo "Downloading GG2 sequences for Naive Bayes classifier..."
wget -N --content-disposition 'https://zenodo.org/records/14169078/files/gg2_2024_09_toSpecies_trainset.fa.gz?download=1'

echo "Downloading SILVA sequences for Naive Bayes classifier..."
wget -N --content-disposition 'https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toSpecies_trainset.fa.gz?download=1'

echo "Downloading GTDB sequences for Naive Bayes classifier..."
wget -N --content-disposition 'https://zenodo.org/records/13984843/files/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz?download=1'

echo "Downloading SILVA sequences and taxonomies for VSEARCH..."
wget -N --content-disposition 'https://data.qiime2.org/2024.10/common/silva-138-99-seqs.qza'
wget -N --content-disposition 'https://data.qiime2.org/2024.10/common/silva-138-99-tax.qza'

echo "Downloading GTDB database for VSEARCH..."
wget -N --content-disposition --no-check-certificate 'https://zenodo.org/record/6912512/files/GTDB_ssu_all_r207.taxonomy.qza?download=1'
wget -N --content-disposition --no-check-certificate 'https://zenodo.org/record/6912512/files/GTDB_ssu_all_r207.qza?download=1'

echo "Downloading GreenGenes2 database for VSEARCH..."
wget -N --content-disposition 'http://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.full-length.fna.qza'
wget -N --content-disposition 'http://ftp.microbio.me/greengenes_release/2024.09/2024.09.backbone.tax.qza'


# Below commands are used to create QZA file for use with VSEARCH with QIIME2. The raw TSV files are downloaded from Zenodo provided by DADA2 for GG2 and GTDB respectively
# GTDB
wget 'https://data.gtdb.ecogenomic.org/releases/release220/220.0/genomic_files_all/ssu_all_r220.fna.gz'
zgrep "^>" ssu_all_r220.fna.gz |  awk '{print gensub(/^>(.*)/, "\\1", "g", $1),gensub(/^>.* (d__.*) \[.*\[.*\[.*/, "\\1", "g", $0)}' OFS=$'\t' > ssu_all_r220.taxonomy.tsv
# Activate QIIME 2 to import
conda activate $HOME/nf_conda/qiime2-amplicon-2024.10-py310-ubuntu-conda-4f67525fc4d0f8938bb2a83db2477cf8
qiime tools import --type 'FeatureData[Sequence]' --input-path ssu_all_r220.fna --output-path GTDB_ssu_all_r220.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path ssu_all_r220.taxonomy.tsv --output-path GTDB_ssu_all_r220.taxonomy.qza

# For GG2
./gg2_to_tax.sh databases/gg2_2024_09_toSpecies_trainset.fa
qiime tools import --type 'FeatureData[Sequence]' --input-path processed.fasta --output-path gg2_2024_09.qza
qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat --input-path output.tsv --output-path gg2_2024_09.taxonomy.qza