nextflow.enable.dsl = 2

include { download_gtdb_db } from '../modules/utils'
include { download_silva_db } from '../modules/utils'
include { download_gg2_db } from '../modules/utils'

include { download_eukaryome_db as download_eukaryome_its_db } from '../modules/utils'
include { download_eukaryome_db as download_eukaryome_18s_db } from '../modules/utils'
include { download_eukaryome_db as download_eukaryome_28s_db } from '../modules/utils'

workflow DOWNLOAD_DATABASES {
    take:
    db_manifest
    requested_dbs

    main:

    if (requested_dbs.contains('gtdb')) {
        download_gtdb_db(
            db_manifest.gtdb.nb.url,
            db_manifest.gtdb.nb.filename,
            db_manifest.gtdb.vsearch.seq_url
        )
    }

    if (requested_dbs.contains('silva')) {
        download_silva_db(
            db_manifest.silva.nb.url,
            db_manifest.silva.nb.filename,
            db_manifest.silva.vsearch.seq_url
        )
    }

    if (requested_dbs.contains('gg2')) {
        download_gg2_db(
            db_manifest.gg2.nb.url,
            db_manifest.gg2.nb.filename,
            db_manifest.gg2.vsearch.seq_url,
            db_manifest.gg2.vsearch.tax_url
        )
    }

    if (requested_dbs.contains('eukaryome_its')) {
        download_eukaryome_its_db(
            'eukaryome_its',
            db_manifest.eukaryome_its.nb.url,
            db_manifest.eukaryome_its.nb.filename,
            db_manifest.eukaryome_its.vsearch.seq_url
        )
    }

    if (requested_dbs.contains('eukaryome_18s')) {
        download_eukaryome_18s_db(
            'eukaryome_18s',
            db_manifest.eukaryome_18s.nb.url,
            db_manifest.eukaryome_18s.nb.filename,
            db_manifest.eukaryome_18s.vsearch.seq_url
        )
    }

    if (requested_dbs.contains('eukaryome_28s')) {
        download_eukaryome_28s_db(
            'eukaryome_28s',
            db_manifest.eukaryome_28s.nb.url,
            db_manifest.eukaryome_28s.nb.filename,
            db_manifest.eukaryome_28s.vsearch.seq_url
        )
    }
}
