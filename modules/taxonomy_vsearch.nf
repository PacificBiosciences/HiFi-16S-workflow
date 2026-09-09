process taxonomy_vsearch_assign {
    label 'highcpu'

    conda (params.enable_conda ? "$projectDir/env/vsearch.yml" : null)
    container "quay.io/biocontainers/vsearch:2.30.0--hd6d6fdc_0"

    publishDir "${params.outdir}/vsearch_tax",
        mode: params.publish_dir_mode

    tag "${db_name}"

    input:
    tuple path(asv_fasta),
          val(db_name),
          path(vsearch_fasta)

    output:
    tuple val(db_name),
          path("${db_name}.vsearch_lca.tsv"),
          emit: vsearch_tax

    tuple val(db_name),
          path("${db_name}.vsearch_hits.tsv"),
          emit: vsearch_hits

    tuple val(db_name),
          path("${db_name}.vsearch_blast6.tsv"),
          emit: vsearch_blast6

    script:
    """
    vsearch \
        --usearch_global "${asv_fasta}" \
        --db "${vsearch_fasta}" \
        --gzip_decompress \
        --id ${params.vsearch_lca_id} \
        --maxaccepts ${params.vsearch_lca_maxaccepts} \
        --maxrejects ${params.vsearch_lca_maxrejects} \
        --query_cov ${params.vsearch_lca_query_cov} \
        --top_hits_only \
        --output_no_hits \
        --strand both \
        --lca_cutoff ${params.vsearch_lca_cutoff} \
        --n_mismatch \
        --notrunclabels \
        --threads ${task.cpus} \
        --userout "${db_name}.vsearch_hits.tsv" \
        --userfields query+target+id+qcov+alnlen+ql+tl \
        --blast6out "${db_name}.vsearch_blast6.tsv" \
        --lcaout "${db_name}.vsearch_lca.tsv"
    """
}


process taxonomy_vsearch_merge {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    input:
    path taxonomy
    path hits
    path asv_table

    output:
    path "vsearch_tax_merged_no_id.tsv", emit: merged_no_id

    script:
    """
    merge_vsearch_taxonomy_with_table.R \\
      ${taxonomy} \\
      ${hits} \\
      ${asv_table} \\
      vsearch_tax_merged_no_id.tsv
    """
}
