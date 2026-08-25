process taxonomy_nb_assign {
    label 'highcpu'

    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/nb_tax", mode: params.publish_dir_mode

    input:
    tuple path(asv_fasta), val(db_name), path(db_fasta)

    output:
    tuple val(db_name), path("${db_name}_nb.tsv"), emit: nb_tax

    script:
    """
    taxonomy_nb_assign.R \\
      ${asv_fasta} \\
      ${db_fasta} \\
      ${db_name} \\
      ${task.cpus}

    harmonise_taxonomy.R \\
      ${db_name}_nb.tsv \\
      ${db_name} \\
      ${db_name}_nb.harmonised.tmp.tsv

    mv ${db_name}_nb.harmonised.tmp.tsv ${db_name}_nb.tsv
    """
}


process taxonomy_nb_best {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    publishDir "${params.outdir}/nb_tax", mode: params.publish_dir_mode

    input:
    path nb_tax_files
    val db_priority

    output:
    path "best_nb_taxonomy.tsv", emit: best_nb_tax
    path "best_nb_taxonomy_withDB.tsv", emit: best_nb_tax_with_db

    script:
    """
    taxonomy_best.R \\
      ${nb_tax_files.join(' ')} \\
      ${db_priority}

    mv best_taxonomy.tsv best_nb_taxonomy.tsv
    mv best_taxonomy_withDB.tsv best_nb_taxonomy_withDB.tsv
    """
}


process taxonomy_nb_merge {
    conda (params.enable_conda ? "$projectDir/env/Rdata_table.yml" : null)
    container "quay.io/biocontainers/r-data.table:1.12.2"

    input:
    path taxonomy
    path asv_table

    output:
    path "best_nb_tax_merged_no_id.tsv", emit: merged_no_id

    script:
    """
    merge_taxonomy_with_table.R \\
      ${taxonomy} \\
      ${asv_table} \\
      best_nb_tax_merged_no_id.tsv
    """
}


