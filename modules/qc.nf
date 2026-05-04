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

process QC_fastq {
    conda (params.enable_conda ? "$projectDir/env/qc_fastq.yml" : null)
    container "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    label 'cpu8'

    publishDir "${params.outdir}/filtered_input_FASTQ", pattern: '*filterQ*.fastq.gz', mode: params.publish_dir_mode
    publishDir "${params.outdir}/qc_raw", pattern: '*.tsv', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path "${sampleID}.seqkit.readstats.tsv", emit: all_seqkit_stats
    path "${sampleID}.seqkit.summarystats.tsv", emit: all_seqkit_summary
    tuple val(sampleID), path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq
    path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq_files

    script:
    """
    seqkit fx2tab -j ${task.cpus} -q --gc -l -H -n -i ${sampleFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \\
        > ${sampleID}.seqkit.readstats.tsv

    seqkit stats -T -j ${task.cpus} -a ${sampleFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \\
        > ${sampleID}.seqkit.summarystats.tsv

    seqkit seq -j ${task.cpus} --min-qual ${params.filterQ} ${sampleFASTQ} \\
        --out-file ${sampleID}.filterQ${params.filterQ}.tmp.fastq.gz

    n_reads=\$(seqkit stats -T ${sampleID}.filterQ${params.filterQ}.tmp.fastq.gz | awk 'NR==2{print \$4}')

    if [ ${params.downsample} -gt 0 ] && [ "\$n_reads" -gt ${params.downsample} ]; then
        seqkit head -n ${params.downsample} \\
            ${sampleID}.filterQ${params.filterQ}.tmp.fastq.gz \\
            --out-file ${sampleID}.filterQ${params.filterQ}.fastq.gz
    else
        mv ${sampleID}.filterQ${params.filterQ}.tmp.fastq.gz \\
          ${sampleID}.filterQ${params.filterQ}.fastq.gz
    fi

    final_reads=\$(seqkit stats -T ${sampleID}.filterQ${params.filterQ}.fastq.gz | awk 'NR==2{print \$4}')

    if [ "\$final_reads" -eq 0 ]; then
        echo "ERROR: QC/downsampling produced 0 reads for sample ${sampleID}" >&2
        echo "Reads after quality filter: \$n_reads" >&2
        echo "Requested downsample: ${params.downsample}" >&2
        exit 1
    fi
    """
}

process cutadapt {
    conda (params.enable_conda ? "$projectDir/env/cutadapt.yml" : null)
    container "quay.io/biocontainers/cutadapt:5.2--py313h8c92656_1"
    label 'cpu16'

    publishDir "${params.outdir}/trimmed_primers_FASTQ", pattern: '*.fastq.gz', mode: params.publish_dir_mode
    publishDir "${params.outdir}/cutadapt_report", pattern: '*.cutadapt.report', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)
    val front_p
    val adapter_p

    output:
    tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq
    tuple val(sampleID), path("${sampleID}.cutadapt.report"), emit: cutadapt_report
    path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq_files

    script:
    """
    cutadapt \
        -g "${front_p}...${adapter_p}" \
        ${sampleFASTQ} \
        -o ${sampleID}.trimmed.fastq.gz \
        -j ${task.cpus} \
        --trimmed-only \
        --revcomp \
        -e 0.1 \
        --json ${sampleID}.cutadapt.report
    """
}

process summarize_cutadapt {
    conda (params.enable_conda ? "$projectDir/env/summarize_cutadapt.yml" : null)
    container "makrezdocker/alpine-jq:1.0"
    label 'cpu_def'

    publishDir "${params.outdir}/cutadapt_summary", pattern: '*.tsv', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(cutadapt_report)

    output:
    path "cutadapt_summary_${sampleID}.tsv", emit: summary_tocollect

    script:
    """
    input_read=\$(jq -r '.read_counts.input' ${cutadapt_report})
    demux_read=\$(jq -r '.read_counts.output' ${cutadapt_report})

    echo -e "sample\tinput_reads\tdemuxed_reads" > cutadapt_summary_${sampleID}.tsv
    echo -e "${sampleID}\t\${input_read}\t\${demux_read}" >> cutadapt_summary_${sampleID}.tsv
    """
}


process QC_fastq_post_trim {
    conda (params.enable_conda ? "$projectDir/env/qc_fastq_post_trim.yml" : null)
    container "quay.io/biocontainers/seqkit:2.13.0--he881be0_0"
    label 'cpu8'

    publishDir "${params.outdir}/qc_post_trim", pattern: '*.post_trim.tsv', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(sampleFASTQ)

    output:
    path "${sampleID}.seqkit.readstats.post_trim.tsv", emit: all_seqkit_stats

    script:
    """
    seqkit fx2tab -j ${task.cpus} -q --gc -l -H -n -i ${sampleFASTQ} | \\
        awk -v sample="${sampleID}" 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"sample"} NR>1{print \$0,sample}' \
        > ${sampleID}.seqkit.readstats.post_trim.tsv
    """
}

process collect_QC {
    conda (params.enable_conda ? "$projectDir/env/collect_qc.yml" : null)
    container "quay.io/biocontainers/csvtk:0.31.0--h9ee0642_0"
    label 'cpu8'

    publishDir "${params.outdir}/reads_QC", mode: params.publish_dir_mode

    input:
    path raw_readstats
    path summarystats
    path cutadapt_stats
    path post_trim_readstats

    output:
    path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
    path "all_samples_seqkit.readstats.post_trim.tsv", emit: all_samples_readstats_post_trim
    path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
    path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
    path "seqkit.summarised_stats.group_by_samples.pretty.tsv", emit: summarised_sample_readstats_pretty
    path "all_samples_cutadapt_stats.tsv", emit: cutadapt_summary

    script:
    """
    csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
    csvtk concat -t -C '%' *.seqkit.readstats.post_trim.tsv > all_samples_seqkit.readstats.post_trim.tsv
    csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv
    csvtk concat -t cutadapt_summary*.tsv > all_samples_cutadapt_stats.tsv

    csvtk summary \
        -t -C '%' \
        -g sample \
        -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median \
        all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv

    csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv \
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
