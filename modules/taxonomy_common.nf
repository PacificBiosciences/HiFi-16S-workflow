process taxonomy_add_md5 {
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)
    container "makrezdocker/alpine-jq:1.0"

    publishDir "${params.outdir}/final", mode: params.publish_dir_mode

    input:
    path tax_table
    val method

    output:
    path "${method}_tax_merged_freq_tax.tsv", emit: tax_with_md5

    script:
    """
    set -euo pipefail

    awk 'BEGIN{FS=OFS="\\t"}
    NR==1{
      print "id", \$0
      next
    }
    {
      cmd = "printf \\"%s\\" \\"" \$1 "\\" | md5sum"
      cmd | getline hashline
      close(cmd)
      split(hashline, a, " ")
      print a[1], \$0
    }' ${tax_table} > ${method}_tax_merged_freq_tax.tsv
    """
}
