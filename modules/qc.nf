/*
===============================================================================

Native QC module for PacBio HiFi 16S preprocessing.

Included:
- metadata inspection
- raw FASTQ QC and quality filtering
- optional cutadapt trimming
- post-trim QC
- QC aggregation

===============================================================================
*/

process inspect_metadata {
    conda (params.enable_conda ? "$projectDir/env/inspect_metadata.yml" : null)
    container "${ workflow.containerEngine == 'singularity' ? 'docker://quay.io/biocontainers/csvtk:0.28.0--h9ee0642_0' : 'quay.io/biocontainers/csvtk:0.28.0--h9ee0642_0' }"
    label 'cpu_def'

    input:
    path sample_sheet
    path metadata

    script:
    """
    set -euo pipefail

    [ -s "${sample_sheet}" ] || { echo "Error: input sample sheet is missing or empty" >&2; exit 1; }
    [ -s "${metadata}" ] || { echo "Error: metadata file is missing or empty" >&2; exit 1; }

    csvtk headers -t "${sample_sheet}" | grep -qx 'sample-id' || { echo "Error: input sample sheet must contain column 'sample-id'" >&2; exit 1; }
    csvtk headers -t "${sample_sheet}" | grep -qx 'filepath' || { echo "Error: input sample sheet must contain column 'filepath'" >&2; exit 1; }
    csvtk headers -t "${metadata}" | grep -qx 'sample_name' || { echo "Error: metadata must contain column 'sample_name'" >&2; exit 1; }

    csvtk cut -t -f sample-id "${sample_sheet}" | tail -n +2 | sort > input.samples.txt
    csvtk cut -t -f sample_name "${metadata}" | tail -n +2 | sort > metadata.samples.txt

    diff input.samples.txt metadata.samples.txt > /dev/null || {
        echo "Error: samples in input sample sheet and metadata do not match" >&2
        exit 1
    }
    """
}

process QC_raw_stats {
    conda (params.enable_conda ? "$projectDir/env/qc_fastq.yml" : null)
    container "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    label 'cpu8'

    publishDir "${params.outdir}/qc_raw",
        pattern: '*.tsv',
        mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path "${sampleID}.seqkit.readstats.tsv", emit: readstats
    path "${sampleID}.seqkit.summarystats.tsv", emit: summarystats

    script:
    """
    set -euo pipefail

    seqkit fx2tab -j ${task.cpus} -q --gc -l -H -n ${sampleFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \\
        > ${sampleID}.seqkit.readstats.tsv

    seqkit stats -T -j ${task.cpus} -a ${sampleFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \\
        > ${sampleID}.seqkit.summarystats.tsv
    """
}


process filter_fastq {
    conda (params.enable_conda ? "$projectDir/env/qc_fastq.yml" : null)
    container "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    label 'cpu8'

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    tuple val(sampleID), path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq

    script:
    """
    set -euo pipefail

    seqkit seq \\
        -j ${task.cpus} \\
        --min-qual ${params.filterQ} \\
        ${sampleFASTQ} \\
        --out-file ${sampleID}.filterQ${params.filterQ}.fastq.gz
    """
}

process downsample_fastq {
    conda (params.enable_conda ? "$projectDir/env/qc_fastq.yml" : null)
    container "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    label 'cpu8'

    publishDir "${params.outdir}/filtered_input_FASTQ",
        pattern: '*.fastq.gz',
        mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(filteredFASTQ)

    output:
    tuple val(sampleID), path("${sampleID}.filterQ${params.filterQ}.downsampled.fastq.gz"), emit: downsampled_fastq
    path "${sampleID}.filterQ${params.filterQ}.downsampled.fastq.gz", emit: downsampled_fastq_files

    script:
    """
    set -euo pipefail

    if [ ${params.downsample} -gt 0 ]; then
        seqkit head \\
            -n ${params.downsample} \\
            ${filteredFASTQ} \\
            --out-file ${sampleID}.filterQ${params.filterQ}.downsampled.fastq.gz
    else
        cp ${filteredFASTQ} ${sampleID}.filterQ${params.filterQ}.downsampled.fastq.gz
    fi
    """
}

process cutadapt {
    conda (params.enable_conda ? "$projectDir/env/cutadapt.yml" : null)
    container "quay.io/biocontainers/cutadapt:5.2--py313h8c92656_1"
    label 'cpu16'

    publishDir "${params.outdir}/cutadapt",
        pattern: '*.log',
        mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(filteredFASTQ)
    val front_p
    val adapter_p

    output:
    tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq
    path "${sampleID}.cutadapt.log", emit: cutadapt_report

    script:
    """
    set -euo pipefail

    cutadapt \\
        -j ${task.cpus} \\
        -g "${front_p}" \\
        -a "${adapter_p}" \\
        -o ${sampleID}.trimmed.fastq.gz \\
        ${filteredFASTQ} \\
        > ${sampleID}.cutadapt.log
    """
}


process summarize_cutadapt {
    conda (params.enable_conda ? "$projectDir/env/summarize_cutadapt.yml" : null)
    container "makrezdocker/alpine-jq:1.0"
    label 'cpu_def'
    label 'cpu_def'

    input:
    path cutadapt_reports

    output:
    path "cutadapt_summary.tsv", emit: summary_tocollect

    script:
    """
    set -euo pipefail

    echo -e "sample\\tcutadapt_log" > cutadapt_summary.tsv

    for log in *.cutadapt.log; do
        sample=\${log%.cutadapt.log}
        echo -e "\${sample}\\t\${log}" >> cutadapt_summary.tsv
    done
    """
}


process QC_fastq_post_trim {
    conda (params.enable_conda ? "$projectDir/env/qc_fastq.yml" : null)
    container "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    label 'cpu8'

    publishDir "${params.outdir}/qc_post_trim",
        pattern: '*.tsv',
        mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(trimmedFASTQ)

    output:
    path "${sampleID}.post_trim.seqkit.readstats.tsv", emit: readstats
    path "${sampleID}.post_trim.seqkit.summarystats.tsv", emit: summarystats
    path "${sampleID}.post_trim.seqkit.readstats.tsv", emit: all_seqkit_stats

    script:
    """
    set -euo pipefail

    seqkit fx2tab -j ${task.cpus} -q --gc -l -H -n ${trimmedFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \\
        > ${sampleID}.post_trim.seqkit.readstats.tsv

    seqkit stats -T -j ${task.cpus} -a ${trimmedFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \\
        > ${sampleID}.post_trim.seqkit.summarystats.tsv
    """
}


process collect_QC {
    conda (params.enable_conda ? "$projectDir/env/collect_qc.yml" : null)
    container "quay.io/biocontainers/csvtk:0.31.0--h9ee0642_0"
    label 'cpu8'

    publishDir "${params.outdir}/reads_QC",
        mode: params.publish_dir_mode

    input:
    path raw_readstats
    path raw_summarystats
    path cutadapt_summary
    path post_trim_readstats

    output:
    path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
    path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
    path "cutadapt_summary.tsv", emit: cutadapt_summary
    path "all_samples_post_trim_seqkit.readstats.tsv", emit: all_samples_post_trim_readstats
    path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
    path "seqkit.summarised_stats.group_by_samples.pretty.tsv", emit: summarised_sample_readstats_pretty

    script:
    """
    set -euo pipefail

    csvtk concat -t -C '%' *.seqkit.readstats.tsv \\
        > all_samples_seqkit.readstats.tsv

    csvtk concat -t -C '%' *.seqkit.summarystats.tsv \\
        > all_samples_seqkit.summarystats.tsv

    cp ${cutadapt_summary} cutadapt_summary.tsv

    csvtk concat -t -C '%' *.post_trim.seqkit.readstats.tsv \\
        > all_samples_post_trim_seqkit.readstats.tsv

    csvtk summary \\
        -t -C '%' \\
        -g sample \\
        -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median \\
        all_samples_post_trim_seqkit.readstats.tsv \\
        > seqkit.summarised_stats.group_by_samples.tsv

    csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv \\
        > seqkit.summarised_stats.group_by_samples.pretty.tsv
    """
}

process collect_QC_skip_cutadapt {
    conda (params.enable_conda ? "$projectDir/env/collect_qc.yml" : null)
    container "quay.io/biocontainers/csvtk:0.31.0--h9ee0642_0"
    label 'cpu8'

    publishDir "${params.outdir}/reads_QC", mode: params.publish_dir_mode

    input:
    path raw_readstats
    path summarystats

    output:
    path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
    path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
    path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
    path "seqkit.summarised_stats.group_by_samples.pretty.tsv", emit: summarised_sample_readstats_pretty

    script:
    """
    csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
    csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv

    csvtk summary \
        -t -C '%' \
        -g sample \
        -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median \
        all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv

    csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv \
        > seqkit.summarised_stats.group_by_samples.pretty.tsv
    """
}
