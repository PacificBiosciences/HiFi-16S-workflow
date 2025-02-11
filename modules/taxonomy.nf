process class_tax {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    cpus params.vsearch_cpu

    input:
    path asv_seq
    path asv_freq
    path vsearch_db
    path vsearch_tax

    output:
    path "taxonomy.vsearch.qza", emit: tax_vsearch
    path "tax_export/taxonomy.tsv", emit: tax_tsv
    path "merged_freq_tax.qzv", emit: tax_freq_tab
    path "vsearch_merged_freq_tax.tsv", emit: tax_freq_tab_tsv
    path "taxonomy.vsearch.blast_results.qza", emit: tax_vsearch_blast

    script:
    """
    qiime feature-classifier classify-consensus-vsearch --i-query $asv_seq \
        --o-classification taxonomy.vsearch.qza \
        --i-reference-reads $vsearch_db \
        --i-reference-taxonomy $vsearch_tax \
        --p-threads $task.cpus \
        --p-maxrejects $params.maxreject \
        --p-maxaccepts $params.maxaccept \
        --p-perc-identity $params.vsearch_identity \
        --p-top-hits-only \
        --o-search-results taxonomy.vsearch.blast_results.qza

    qiime tools export --input-path taxonomy.vsearch.qza --output-path tax_export

    qiime feature-table transpose --i-table $asv_freq \
        --o-transposed-feature-table transposed-asv.qza

    qiime metadata tabulate --m-input-file $asv_seq \
        --m-input-file taxonomy.vsearch.qza \
        --m-input-file transposed-asv.qza \
        --o-visualization merged_freq_tax.qzv

    qiime tools export --input-path merged_freq_tax.qzv \
        --output-path merged_freq_tax_tsv

    mv merged_freq_tax_tsv/metadata.tsv vsearch_merged_freq_tax.tsv
    """
}

process dada2_assignTax {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", pattern: 'best_tax*', mode: params.publish_dir_mode
    publishDir "$params.outdir/nb_tax", mode: params.publish_dir_mode
    cpus params.vsearch_cpu

    input:
    path asv_seq_fasta
    path asv_seq
    path asv_freq
    path silva_db
    path gtdb_db
    path refseq_db
    path assign_script
    val(db_to_prioritize)

    output:
    path "best_tax.qza", emit:best_nb_tax_qza
    path "best_taxonomy.tsv", emit: best_nb_tax
    path "best_taxonomy_withDB.tsv"
    path "best_tax_merged_freq_tax.tsv", emit: best_nb_tax_tsv
    path "silva_nb.tsv"
    path "gtdb_nb.tsv"
    path "gg2_nb.tsv"

    script:
    """
    Rscript --vanilla $assign_script $asv_seq_fasta $task.cpus $silva_db $gtdb_db $refseq_db 80 $db_to_prioritize

    qiime feature-table transpose --i-table $asv_freq \
        --o-transposed-feature-table transposed-asv.qza

    qiime tools import --type "FeatureData[Taxonomy]" \
        --input-format "TSVTaxonomyFormat" \
        --input-path best_taxonomy.tsv --output-path best_tax.qza

    qiime metadata tabulate --m-input-file $asv_seq \
        --m-input-file best_tax.qza \
        --m-input-file transposed-asv.qza \
        --o-visualization merged_freq_tax.qzv

    qiime tools export --input-path merged_freq_tax.qzv \
        --output-path merged_freq_tax_tsv

    mv merged_freq_tax_tsv/metadata.tsv best_tax_merged_freq_tax.tsv
    """
}

process export_biom {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_freq
    path tax_tsv
    path tax_tsv_vsearch

    output:
    path "feature-table-tax.biom", emit:biom
    path "feature-table-tax_vsearch.biom", emit:biom_vsearch

    script:
    """
    qiime tools export --input-path $asv_freq --output-path asv_freq/

    sed 's/Feature ID/#OTUID/' $tax_tsv | sed 's/Taxon/taxonomy/' | \
        sed 's/Consensus/confidence/' > biom-taxonomy.tsv

    biom add-metadata -i asv_freq/feature-table.biom \
        -o feature-table-tax.biom \
        --observation-metadata-fp biom-taxonomy.tsv \
        --sc-separated taxonomy

    sed 's/Feature ID/#OTUID/' $tax_tsv_vsearch | sed 's/Taxon/taxonomy/' | \
        sed 's/Consensus/confidence/' > biom-taxonomy_vsearch.tsv

    biom add-metadata -i asv_freq/feature-table.biom \
        -o feature-table-tax_vsearch.biom \
        --observation-metadata-fp biom-taxonomy_vsearch.tsv \
        --sc-separated taxonomy
    """
}

process export_biom_skip_nb {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/qiime2/amplicon@sha256:4038fd785bf4e76ddd6ec7a7f57abe94cdca6c5cd0a93d0924971a74eabd7cf2"
    publishDir "$params.outdir/results", mode: params.publish_dir_mode
    label 'cpu_def'

    input:
    path asv_freq
    path tax_tsv_vsearch

    output:
    path "feature-table-tax_vsearch.biom", emit:biom_vsearch

    script:
    """
    qiime tools export --input-path $asv_freq --output-path asv_freq/

    sed 's/Feature ID/#OTUID/' $tax_tsv_vsearch | sed 's/Taxon/taxonomy/' | \
        sed 's/Consensus/confidence/' > biom-taxonomy_vsearch.tsv

    biom add-metadata -i asv_freq/feature-table.biom \
        -o feature-table-tax_vsearch.biom \
        --observation-metadata-fp biom-taxonomy_vsearch.tsv \
        --sc-separated taxonomy
    """
}