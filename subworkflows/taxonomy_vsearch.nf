nextflow.enable.dsl = 2

include {
    taxonomy_vsearch_assign
    taxonomy_vsearch_merge
} from '../modules/taxonomy_vsearch'

include {
    taxonomy_add_md5 as taxonomy_vsearch_add_md5
} from '../modules/taxonomy_common'


workflow TAXONOMY_VSEARCH {

    take:
    vsearch_inputs_ch
    asv_table_tsv

    main:

    taxonomy_vsearch_assign(
        vsearch_inputs_ch
    )

    taxonomy_vsearch_merge(
        taxonomy_vsearch_assign.out.vsearch_tax
            .map { db_name, taxonomy -> taxonomy },

        taxonomy_vsearch_assign.out.vsearch_hits
            .map { db_name, hits -> hits },

        asv_table_tsv
    )

    taxonomy_vsearch_add_md5(
        taxonomy_vsearch_merge.out.merged_no_id,
        'vsearch'
    )

    emit:
    vsearch_tax         = taxonomy_vsearch_assign.out.vsearch_tax
    vsearch_hits        = taxonomy_vsearch_assign.out.vsearch_hits
    final_vsearch_table = taxonomy_vsearch_add_md5.out.tax_with_md5
}
