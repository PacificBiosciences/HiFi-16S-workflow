process taxonomy_vsearch_assign {
    label 'highcpu'

    conda (params.enable_conda ? "$projectDir/env/vsearch.yml" : null)
    container "quay.io/biocontainers/vsearch:2.30.0--hd6d6fdc_0"

    publishDir "${params.outdir}/vsearch_tax", mode: params.publish_dir_mode

    tag "${db_name}"

    input:
    tuple path(asv_fasta),
          val(db_name),
          path(vsearch_fasta)

    output:
    tuple val(db_name),
          path("${db_name}.vsearch_lca.tsv"),
          emit: vsearch_tax

    script:
    """
    vsearch \
        --usearch_global "${asv_fasta}" \
        --db "${vsearch_fasta}" \
        --gzip_decompress \
        --id ${params.vsearch_lca_id} \
        --query_cov ${params.vsearch_lca_query_cov} \
        --maxaccepts ${params.vsearch_lca_maxaccepts} \
        --maxrejects ${params.vsearch_lca_maxrejects} \
        --top_hits_only \
        --output_no_hits \
        --lca_cutoff ${params.vsearch_lca_cutoff} \
        --n_mismatch \
        --notrunclabels \
        --threads ${task.cpus} \
        --lcaout "${db_name}.vsearch_lca.tsv"
    """
}
