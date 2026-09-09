process download_nb_db {

    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir { "${params.db_base_dir}/${db_name}/nb" },
        pattern: "*.fa.gz",
        mode: "copy"

    label 'cpu_def'

    input:
    val db_name
    val nb_url
    val nb_filename

    output:
    tuple val(db_name),
          path("${nb_filename}"),
          emit: nb

    script:
    """
    set -euo pipefail

    echo "Downloading ${db_name} Naive Bayes database..."

    wget \\
        -O ${nb_filename} \\
        "${nb_url}"
    """
}

process download_eukaryome_db {

    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.2"

    publishDir { "${params.db_base_dir}/${db_name}" },
        mode: "copy"

    label 'cpu_def'

    input:
    val db_name
    val nb_url
    val nb_filename
    val vsearch_url
    val vsearch_filename

    output:
    tuple val(db_name),
          path("nb/${nb_filename}"),
          emit: nb

    tuple val(db_name),
          path("vsearch/${vsearch_filename}"),
          emit: vsearch

    script:
    """
    mkdir -p nb vsearch

    prepare_eukaryome_db.sh \
        "${db_name}" \
        "${nb_url}" \
        "nb/${nb_filename}" \
        "${vsearch_url}" \
        "vsearch/${vsearch_filename}"
    """
}
