nextflow.enable.dsl = 2


include {
    taxonomy_vsearch_assign
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


    emit:
    vsearch_tax         = taxonomy_vsearch_assign.out.vsearch_tax
}
