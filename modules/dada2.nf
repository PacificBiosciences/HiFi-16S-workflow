process dada2_filter_ccs {
    
	label 'cpu_def'

    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/filter_stats", pattern: '*.filter_stats.tsv', mode: params.publish_dir_mode

    input:
    tuple val(sampleID), path(trimmed_fastq)

    output:
    tuple val(sampleID), path("${sampleID}.filtered.fastq.gz"), emit: filtered_fastq
    path "${sampleID}.filter_stats.tsv", emit: filter_stats

    script:
    def minLenArg = params.min_len != null ? params.min_len : 'NA'
    def maxLenArg = params.max_len != null ? params.max_len : 'NA'
    def maxEeArg  = params.max_ee  != null ? params.max_ee  : 'NA'

    """
    filter_and_trim.R \\
      ${trimmed_fastq} \\
      ${sampleID}.filtered.fastq.gz \\
      ${sampleID}.filter_stats.tsv \\
      ${maxEeArg} \\
      ${minLenArg} \\
      ${maxLenArg} \\
      0 \\
      FALSE
    """
}

process subsample_for_error_model {
    label 'highparallel'
    conda (params.enable_conda ? "$projectDir/env/qc_fastq.yml" : null)
    container "quay.io/biocontainers/seqtk:1.4--he4a0461_2"

    input:
    tuple val(sampleID), path(filtered_fastq)

    output:
    path "${sampleID}.error_sample.fastq.gz", emit: sampled_fastq

    script:
    """
    seqtk sample \
      ${filtered_fastq} \
      ${params.error_model_reads_per_sample} \
      | gzip -c > ${sampleID}.error_sample.fastq.gz
    """
}

process concatenate_error_model_reads {
    container "makrezdocker/alpine-jq:1.0"
    conda (params.enable_conda ? "$projectDir/env/jq.yml" : null)

    label 'lowcpu'

    publishDir "${params.outdir}/dada2/error_model",
        pattern: 'error_model_sample.fastq.gz',
        mode: params.publish_dir_mode

    input:
    path sampled_fastqs

    output:
    path "error_model_sample.fastq.gz", emit: error_sample

    script:
    """
    cat ${sampled_fastqs.join(' ')} > error_model_sample.fastq.gz
    """
}

process learn_errors {
    label 'highcpu'
    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/error_model", pattern: 'errorfun.rds', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/error_model", pattern: 'plot_error_model.pdf', mode: params.publish_dir_mode

    input:
    path error_sample

    output:
    path "errorfun.rds", emit: error_model
    path "plot_error_model.pdf", emit: error_plot

    script:
    """
    learn_errors.R \\
      ${params.learn_nbases ?: '1e9'} \\
      ${params.binned_quality_scores} \\
      ${error_sample}
    """
}

process dada2_denoise_independent {
    label 'highparallel'
    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/denoised", pattern: '*.dada.rds', mode: params.publish_dir_mode


    input:
    tuple val(sampleID), path(filtered_fastq), path(error_model)

    output:
    tuple val(sampleID), path("${sampleID}.dada.rds"), emit: dada_rds
    path "${sampleID}.denoise_stats.tsv", emit: denoise_stats

    script:
    def omegaCArg = params.omegac != null ? params.omegac : '1e-40'
    def bandSizeArg = params.band_size != null ? params.band_size : '16'
    def hpGapArg = params.homopolymer_gap_penalty != null ? params.homopolymer_gap_penalty : 'NA'

    """
    denoise_independent.R \\
      ${filtered_fastq} \\
      ${error_model} \\
      ${sampleID}.dada.rds \\
      ${omegaCArg} \\
      ${bandSizeArg} \\
      ${hpGapArg}
    """
}

process dada2_make_seqtab {
	label 'cpu16'
    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/seqtab", pattern: 'seqtab.rds', mode: params.publish_dir_mode

    input:
    path dada_rds_files

    output:
    path "seqtab.rds", emit: seqtab_rds

    script:
    """
    make_seqtab.R \\
      ${dada_rds_files.join(' ')}
    """
}

process dada2_remove_chimeras {
    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/chimera_removal", pattern: 'seqtab_nochim.rds', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/chimera_removal", pattern: 'chimera_stats.tsv', mode: params.publish_dir_mode

    input:
    path seqtab_rds

    output:
    path "seqtab_nochim.rds", emit: seqtab_nochim_rds
    path "chimera_stats.tsv", emit: chimera_stats

    script:
    def chimeraMethodArg = params.chimera_method != null ? params.chimera_method : 'consensus'
    def minParentFoldArg = params.min_parent_fold != null ? params.min_parent_fold : '1.0'

    """
    remove_chimeras.R \\
      ${seqtab_rds} \\
      seqtab_nochim.rds \\
      chimera_stats.tsv \\
      ${chimeraMethodArg} \\
      ${minParentFoldArg}
    """
}

process dada2_filter_asvs {
    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'seqtab_nochim_filtered.rds', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'dada2_ASV.fasta', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'dada2_table_filtered.tsv', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'filter_asv_stats.tsv', mode: params.publish_dir_mode

    input:
    path seqtab_nochim_rds
    val min_asv_totalfreq
    val min_asv_sample

    output:
    path "seqtab_nochim_filtered.rds", emit: seqtab_filtered_rds
    path "dada2_ASV.fasta", emit: asv_fasta
    path "dada2_table_filtered.tsv", emit: asv_table_tsv
    path "filter_asv_stats.tsv", emit: filter_asv_stats

    script:
    """
    filter_asvs.R \\
      ${seqtab_nochim_rds} \\
      seqtab_nochim_filtered.rds \\
      dada2_table_filtered.tsv \\
      dada2_ASV.fasta \\
      filter_asv_stats.tsv \\
      ${min_asv_totalfreq} \\
      ${min_asv_sample}
    """
}

process dada2_stats {
    conda (params.enable_conda ? "$projectDir/env/dada2.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/final", mode: params.publish_dir_mode

    cpus 1

    input:
    path filter_stats_files
    path denoise_stats_files
    path seqtab_nochim_rds
    path seqtab_filtered_rds
    path metadata

    output:
    path "dada2_tracking.tsv", emit: dada2_tracking_tsv
    path "frequency.tsv", emit: frequency_tsv
    path "rarefaction_depth_suggested.txt", emit: rarefaction_depth_file

    script:
    """
    dada2_stats.R \\
      ${seqtab_nochim_rds} \\
      ${seqtab_filtered_rds} \\
      ${metadata} \\
      dada2_tracking.tsv \\
      frequency.tsv \\
      rarefaction_depth_suggested.txt \\
      ${filter_stats_files.join(' ')} \\
      --DENOISE_STATS-- \\
      ${denoise_stats_files.join(' ')}
    """
}
