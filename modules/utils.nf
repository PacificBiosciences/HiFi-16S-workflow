process download_gtdb_db {
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.db_base_dir}/gtdb/nb", pattern: "*.fa.gz", mode: "copy"
    publishDir "${params.db_base_dir}/gtdb/vsearch", pattern: "sequences.fasta", mode: "copy"
    publishDir "${params.db_base_dir}/gtdb/vsearch", pattern: "taxonomy.tsv", mode: "copy"

    label 'cpu_def'

    input:
    val nb_url
    val nb_filename
    val vsearch_seq_url

    output:
    path "${nb_filename}"
    path "sequences.fasta"
    path "taxonomy.tsv"

    script:
    """
    set -euo pipefail

    wget -O ${nb_filename} "${nb_url}"
    wget -O gtdb_sequences.fna.gz "${vsearch_seq_url}"

    gunzip -c gtdb_sequences.fna.gz > sequences.fasta

    gunzip -c gtdb_sequences.fna.gz | \\
      grep '^>' | \\
      awk 'BEGIN{OFS="\\t"}
      {
        id=\$1
        sub(/^>/,"",id)

        tax=\$0
        sub(/^>[^ ]+[[:space:]]+/, "", tax)
        sub(/ \\[.*/, "", tax)

        print id, tax
      }' > taxonomy.tsv

    rm -f gtdb_sequences.fna.gz
    """
}

process download_silva_db {
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.db_base_dir}/silva/nb", pattern: "*.fa.gz", mode: "copy"
    publishDir "${params.db_base_dir}/silva/vsearch", pattern: "sequences.fasta", mode: "copy"
    publishDir "${params.db_base_dir}/silva/vsearch", pattern: "taxonomy.tsv", mode: "copy"

    label 'cpu_def'

    input:
    val nb_url
    val nb_filename
    val vsearch_seq_url

    output:
    path "${nb_filename}"
    path "sequences.fasta"
    path "taxonomy.tsv"

    script:
    """
    set -euo pipefail

    wget -O ${nb_filename} "${nb_url}"
    wget -O silva_sequences.fasta.gz "${vsearch_seq_url}"

    gunzip -c silva_sequences.fasta.gz > sequences.fasta

    gunzip -c silva_sequences.fasta.gz | \\
      grep '^>' | \\
      awk 'BEGIN{OFS="\\t"}
      {
        id=\$1
        sub(/^>/,"",id)

        tax=\$0
        sub(/^>[^ ]+[[:space:]]+/, "", tax)

        print id, tax
      }' > taxonomy.tsv

    rm -f silva_sequences.fasta.gz
    """
}

process download_gg2_db {
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.db_base_dir}/gg2/nb", pattern: "*.fa.gz", mode: "copy"
    publishDir "${params.db_base_dir}/gg2/vsearch", pattern: "sequences.fasta", mode: "copy"
    publishDir "${params.db_base_dir}/gg2/vsearch", pattern: "taxonomy.tsv", mode: "copy"

    label 'cpu_def'

    input:
    val nb_url
    val nb_filename
    val vsearch_seq_url
    val vsearch_tax_url

    output:
    path "${nb_filename}"
    path "sequences.fasta"
    path "taxonomy.tsv"

    script:
    """
    set -euo pipefail

    echo "Downloading GG2 NB database..."
    wget -O ${nb_filename} "${nb_url}"

    echo "Downloading GG2 VSEARCH sequences..."
    wget -O gg2_sequences.fna.gz "${vsearch_seq_url}"

    echo "Downloading GG2 taxonomy..."
    wget -O gg2_taxonomy.tsv.gz "${vsearch_tax_url}"

    gunzip -c gg2_sequences.fna.gz > sequences.fasta

    echo "Converting GG2 taxonomy to VSEARCH format..."

    gunzip -c gg2_taxonomy.tsv.gz | \\
      awk -F '\\t' 'BEGIN{OFS="\\t"}
      NR==1 && \$1 ~ /Feature ID/ {next}
      {
        print \$1, \$2
      }' > taxonomy.tsv

    rm -f gg2_sequences.fna.gz gg2_taxonomy.tsv.gz
    """
}
