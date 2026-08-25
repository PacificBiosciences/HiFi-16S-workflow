process taxonomy_vsearch_assign {
    label 'highcpu'

    conda (params.enable_conda ? "$projectDir/env/vsearch.yml" : null)
    container "quay.io/biocontainers/vsearch:2.8.0--hfc679d8_0"

    publishDir "${params.outdir}/vsearch_tax", mode: params.publish_dir_mode

    input:
    tuple path(asv_fasta), val(db_name), path(vsearch_fasta), path(vsearch_taxonomy)

    output:
    tuple val(db_name), path("${db_name}_vsearch.tsv"), emit: vsearch_tax

    script:
    """
    taxonomy_vsearch_assign.sh \\
      ${asv_fasta} \\
      ${vsearch_fasta} \\
      ${vsearch_taxonomy} \\
      ${db_name} \\
      ${task.cpus} \\
      ${params.maxreject} \\
      ${params.maxaccept} \\
      ${params.vsearch_identity} \\
      ${db_name}_vsearch.tsv
    """
}

process taxonomy_vsearch_merge {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    publishDir "${params.outdir}/vsearch_tax", mode: params.publish_dir_mode

    input:
    tuple val(db_name), path(taxonomy), path(asv_table)

    output:
    tuple val(db_name), path("${db_name}_vsearch_counts.tsv"), emit: vsearch_tax_counts

    script:
    """
    merge_vsearch_taxonomy_with_table.R \\
      ${taxonomy} \\
      ${asv_table} \\
      ${db_name}_vsearch_counts.tsv
    """
}


process taxonomy_vsearch_best {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    publishDir "${params.outdir}/vsearch_tax", mode: params.publish_dir_mode

    input:
    path vsearch_tax_files
    val db_priority

    output:
    path "best_vsearch_taxonomy.tsv", emit: best_vsearch_tax
    path "best_vsearch_taxonomy_withDB.tsv", emit: best_vsearch_tax_with_db

    script:
    """
    taxonomy_vsearch_best.R \\
      ${vsearch_tax_files.join(' ')} \\
      ${db_priority}
    """
}

process taxonomy_vsearch_consensus {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    input:
    tuple val(db_name), path(vsearch_hits)

    output:
    tuple val(db_name), path("${db_name}.vsearch.consensus.tsv"), emit: vsearch_consensus

    script:
    """
    vsearch_consensus.R \\
      ${vsearch_hits} \\
      ${db_name}.vsearch.consensus.tsv \\
      ${params.vsearch_min_consensus}
    """
}
