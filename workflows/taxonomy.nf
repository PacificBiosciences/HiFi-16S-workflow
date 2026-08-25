nextflow.enable.dsl = 2

include {
    TAXONOMY_NB
} from '../subworkflows/taxonomy_nb'

include {
    TAXONOMY_VSEARCH
} from '../subworkflows/taxonomy_vsearch'


workflow TAXONOMY_WORKFLOW {

    take:
    db_manifest
    asv_fasta
    asv_table_tsv

    main:

    selected_nb_dbs = params.nb_databases instanceof String
        ? params.nb_databases.split(',')*.trim().findAll()
        : params.nb_databases

    nb_db_defs = selected_nb_dbs.collect { db ->
        if (!db_manifest.containsKey(db)) {
            error "Database '${db}' is listed in params.nb_databases but not found in ${params.databases_yaml}"
        }

        if (!db_manifest[db].containsKey('nb')) {
            error "Database '${db}' has no 'nb' section in ${params.databases_yaml}"
        }

        tuple(
            db,
            file(
                "${params.db_base_dir}/${db}/nb/${db_manifest[db].nb.filename}",
                checkIfExists: true
            )
        )
    }

    nb_db_ch = Channel.fromList(nb_db_defs)

    nb_inputs_ch = asv_fasta
        .combine(nb_db_ch)
        .map { fasta, db_name, db_fasta ->
            tuple(fasta, db_name, db_fasta)
        }

    TAXONOMY_NB(
        nb_inputs_ch,
        asv_table_tsv,
        params.db_to_prioritize
    )

    selected_vsearch_dbs = params.vsearch_databases instanceof String
        ? params.vsearch_databases.split(',')*.trim().findAll()
        : params.vsearch_databases

    vsearch_db_defs = selected_vsearch_dbs.collect { db ->
        if (!db_manifest.containsKey(db)) {
            error "Database '${db}' is listed in params.vsearch_databases but not found in ${params.databases_yaml}"
        }

        if (!db_manifest[db].containsKey('vsearch')) {
            error "Database '${db}' has no 'vsearch' section in ${params.databases_yaml}"
        }

        tuple(
            db,
            file(
                "${params.db_base_dir}/${db}/vsearch/${db_manifest[db].vsearch.seq_filename}",
                checkIfExists: true
            ),
            file(
                "${params.db_base_dir}/${db}/vsearch/${db_manifest[db].vsearch.tax_filename}",
                checkIfExists: true
            )
        )
    }

    vsearch_db_ch = Channel.fromList(vsearch_db_defs)

    vsearch_inputs_ch = asv_fasta
        .combine(vsearch_db_ch)
        .map { fasta, db_name, vsearch_fasta, vsearch_taxonomy ->
            tuple(fasta, db_name, vsearch_fasta, vsearch_taxonomy)
        }

    TAXONOMY_VSEARCH(
        vsearch_inputs_ch,
        asv_table_tsv,
        params.db_to_prioritize
    )
}
