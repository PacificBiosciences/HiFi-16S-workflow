def parseRequestedDbs(download_targets) {
    download_targets ? download_targets.tokenize(',')*.trim() : []
}

def validateDownloadParams(params, requested_dbs, db_manifest) {
    if (!params.download_targets) {
        error "When --download_db true, you must also provide --download_targets"
    }

    def valid_dbs = [
        'silva',
        'gtdb',
        'gg2',
        'euk_ssu',
        'euk_lsu',
        'euk_long',
        'euk_its'
    ]

    requested_dbs.each { db ->
        if (!valid_dbs.contains(db)) {
            error "Invalid database '${db}'. Allowed values: ${valid_dbs.join(', ')}"
        }

        if (!db_manifest[db]) {
            error "Database '${db}' is missing from conf/databases.yml"
        }
    }

    if (!params.db_base_dir) {
        error "Missing required parameter: --db_base_dir"
    }
}

def validatePreprocessParams(params) {
    if (!params.input) {
        error "Missing required parameter: --input"
    }
    if (!params.metadata) {
        error "Missing required parameter: --metadata"
    }
    if (!file(params.input).exists()) {
        error "Input sample sheet not found: ${params.input}"
    }
    if (!file(params.metadata).exists()) {
        error "Metadata file not found: ${params.metadata}"
    }

    def n_sample = file(params.input).countLines() - 1
    if (n_sample < 1) {
        error "Input sample sheet appears empty: ${params.input}"
    }

    return n_sample
}

def buildRunLog(params, version, n_sample) {
    def trim_cutadapt = params.skip_primer_trim ? "No" : "Yes"

    return """
Parameters for native pb-16S preprocessing
=========================================
Number of samples in samples TSV: ${n_sample}
Metadata file: ${params.metadata}
Filter input reads above Q: ${params.filterQ}
Downsample reads per sample (0 = disabled): ${params.downsample}
Trim primers with cutadapt: ${trim_cutadapt}
Forward primer: ${params.front_p}
Reverse primer: ${params.adapter_p}
Output directory: ${params.outdir}
Execution profile uses conda: ${params.enable_conda}
Execution profile uses containers: ${params.enable_container}
Pipeline version: ${version}
""".stripIndent()
}
