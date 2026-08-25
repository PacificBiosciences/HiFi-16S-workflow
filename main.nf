nextflow.enable.dsl = 2

include {
    parseRequestedDbs
    validateDownloadParams
    validatePreprocessParams
    buildRunLog
} from './modules/validation'

include { DOWNLOAD_DATABASES } from './workflows/download_databases'
include { PB16S_PREPROCESS }   from './workflows/pb16s_preprocess'
include { DADA2_WORKFLOW }     from './workflows/dada2'
include { TAXONOMY_WORKFLOW }  from './workflows/taxonomy'

def helpMessage() {
    """
===============================================================================
PacBio HiFi 16S preprocessing workflow
===============================================================================

DESCRIPTION
    Native Nextflow implementation of the PacBio HiFi 16S workflow.

USAGE

    Preprocessing:
      nextflow run main.nf \\
        --input sampledata.tsv \\
        --metadata metadata.tsv \\
        --db_base_dir databases \\
        -profile conda

    Download databases:
      nextflow run main.nf \\
        --download_db \\
        --download_targets silva,gg2,gtdb \\
        -profile conda

REQUIRED PARAMETERS

    --input               Sample sheet TSV
    --metadata            Metadata TSV
    --db_base_dir         Database directory

DATABASE DOWNLOAD

    --download_db         Enable database download mode
    --download_targets    Comma-separated list of databases

GENERAL

    --outdir              Output directory
    --filterQ             Minimum CCS read quality
    --downsample          Reads per sample for QC plots
    --error_model         Existing DADA2 error model

PRIMER TRIMMING

    --skip_primer_trim
    --front_p
    --adapter_p

DADA2

    --min_len
    --max_len
    --max_ee
    --learn_nbases
    --band_size
    --homopolymer_gap_penalty
    --omegac

CHIMERA REMOVAL

    --chimera_method
    --min_parent_fold

ASV FILTERING

    --min_asv_total_freq
    --min_asv_sample

TAXONOMY

    --databases_yaml
    --nb_databases
    --vsearch_databases
    --db_to_prioritize

EXECUTION PROFILES

    standard
    conda
    docker
    singularity
    singularity,slurm

HELP

    nextflow run main.nf --help

VERSION

    nextflow run main.nf --version

===============================================================================
"""
}


workflow {

        if (params.help) {
            log.info helpMessage()
            return
        }

        if (params.version) {
            log.info workflow.manifest.version
            return
        }

        def db_manifest_file = file(params.databases_yaml)

        if (!db_manifest_file.exists()) {
            error "Database manifest not found: ${db_manifest_file}"
        }

        def db_manifest = new groovy.yaml.YamlSlurper().parse(db_manifest_file)

        if (params.download_db) {
            def requested_dbs = parseRequestedDbs(params.download_targets)

            validateDownloadParams(params, requested_dbs, db_manifest)

            DOWNLOAD_DATABASES(db_manifest, requested_dbs)
        }
        else {
            def n_sample = validatePreprocessParams(params)

            log.info buildRunLog(
                params,
                workflow.manifest.version,
                n_sample
            )

            sample_sheet_ch = Channel.fromPath(params.input)
            metadata_ch     = Channel.fromPath(params.metadata)

            /*
            * Prepare reads
            */
            PB16S_PREPROCESS(
                sample_sheet_ch,
                metadata_ch
            )

            /*
            * DADA2
            */
            DADA2_WORKFLOW(
                PB16S_PREPROCESS.out.reads_for_dada2,
                metadata_ch
            )

            /*
            * Taxonomy
            */
            TAXONOMY_WORKFLOW(
                db_manifest,
                DADA2_WORKFLOW.out.asv_fasta,
                DADA2_WORKFLOW.out.asv_table_tsv
            )
    }
}
