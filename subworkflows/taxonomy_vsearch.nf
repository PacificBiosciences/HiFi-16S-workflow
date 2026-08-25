nextflow.enable.dsl = 2

include {
    taxonomy_vsearch_assign
    taxonomy_vsearch_consensus
    taxonomy_vsearch_merge
    taxonomy_vsearch_best
} from '../modules/taxonomy_vsearch'

include { taxonomy_add_md5 as taxonomy_vsearch_add_md5 } \
    from '../modules/taxonomy_common'


workflow TAXONOMY_VSEARCH {

    take:
    vsearch_inputs_ch
    asv_table_tsv
    db_priority

    main:

    taxonomy_vsearch_assign(vsearch_inputs_ch)

    taxonomy_vsearch_consensus(
        taxonomy_vsearch_assign.out.vsearch_tax
    )

    vsearch_merge_inputs_ch = taxonomy_vsearch_consensus.out.vsearch_consensus
        .combine(asv_table_tsv)
        .map { db_name, taxonomy, asv_table ->
            tuple(db_name, taxonomy, asv_table)
        }
    taxonomy_vsearch_merge(vsearch_merge_inputs_ch)

    taxonomy_vsearch_best(
        taxonomy_vsearch_merge.out.vsearch_tax_counts
            .map { db_name, table -> table }
            .collect(),
        db_priority
    )

    taxonomy_vsearch_add_md5(
        taxonomy_vsearch_best.out.best_vsearch_tax,
        'best_vsearch'
    )

    emit:
    vsearch_tax              = taxonomy_vsearch_assign.out.vsearch_tax
    vsearch_tax_counts       = taxonomy_vsearch_merge.out.vsearch_tax_counts
    best_vsearch_tax         = taxonomy_vsearch_best.out.best_vsearch_tax
    best_vsearch_tax_with_db = taxonomy_vsearch_best.out.best_vsearch_tax_with_db
    final_vsearch_table      = taxonomy_vsearch_add_md5.out.tax_with_md5
}
