nextflow.enable.dsl = 2


include { download_nb_db as download_gtdb_db } \
    from '../modules/utils'

include { download_nb_db as download_silva_db } \
    from '../modules/utils'

include { download_nb_db as download_gg2_db } \
    from '../modules/utils'


include { download_eukaryome_db as download_euk_ssu_db } \
    from '../modules/utils'

include { download_eukaryome_db as download_euk_lsu_db } \
    from '../modules/utils'

include { download_eukaryome_db as download_euk_long_db } \
    from '../modules/utils'

include { download_eukaryome_db as download_euk_its_db } \
    from '../modules/utils'


workflow DOWNLOAD_DATABASES {

    take:
    db_manifest
    requested_dbs

    main:

    // -------------------------------------------------------------------------
    // 16S databases
    // Naive Bayes only
    // -------------------------------------------------------------------------

    if (requested_dbs.contains('gtdb')) {
        download_gtdb_db(
            'gtdb',
            db_manifest.gtdb.nb.url,
            db_manifest.gtdb.nb.filename
        )
    }


    if (requested_dbs.contains('silva')) {
        download_silva_db(
            'silva',
            db_manifest.silva.nb.url,
            db_manifest.silva.nb.filename
        )
    }


    if (requested_dbs.contains('gg2')) {
        download_gg2_db(
            'gg2',
            db_manifest.gg2.nb.url,
            db_manifest.gg2.nb.filename
        )
    }


    // -------------------------------------------------------------------------
    // EUKARYOME databases
    // Naive Bayes + VSEARCH/SINTAX
    // -------------------------------------------------------------------------

    if (requested_dbs.contains('euk_ssu')) {
        download_euk_ssu_db(
            'euk_ssu',

            db_manifest.euk_ssu.nb.url,
            db_manifest.euk_ssu.nb.filename,

            db_manifest.euk_ssu.vsearch.url,
            db_manifest.euk_ssu.vsearch.filename
        )
    }


    if (requested_dbs.contains('euk_lsu')) {
        download_euk_lsu_db(
            'euk_lsu',

            db_manifest.euk_lsu.nb.url,
            db_manifest.euk_lsu.nb.filename,

            db_manifest.euk_lsu.vsearch.url,
            db_manifest.euk_lsu.vsearch.filename
        )
    }


    if (requested_dbs.contains('euk_long')) {
        download_euk_long_db(
            'euk_long',

            db_manifest.euk_long.nb.url,
            db_manifest.euk_long.nb.filename,

            db_manifest.euk_long.vsearch.url,
            db_manifest.euk_long.vsearch.filename
        )
    }


    if (requested_dbs.contains('euk_its')) {
        download_euk_its_db(
            'euk_its',

            db_manifest.euk_its.nb.url,
            db_manifest.euk_its.nb.filename,

            db_manifest.euk_its.vsearch.url,
            db_manifest.euk_its.vsearch.filename
        )
    }
}
