nextflow.enable.dsl = 2

include {
    dada2_filter_ccs
    subsample_for_error_model
    concatenate_error_model_reads
    learn_errors
    dada2_denoise_independent
    dada2_make_seqtab
    dada2_remove_chimeras
    dada2_filter_asvs
    dada2_stats
} from '../modules/dada2'

workflow DADA2_WORKFLOW {

    take:
    reads_for_dada2
    metadata_ch

    main:

    dada2_filter_ccs(reads_for_dada2)

    if (params.error_model) {
        error_model_ch = Channel.fromPath(params.error_model, checkIfExists: true)
    }
    else {
        subsample_for_error_model(dada2_filter_ccs.out.filtered_fastq)

        concatenate_error_model_reads(
            subsample_for_error_model.out.sampled_fastq.collect()
        )

        learn_errors(
            concatenate_error_model_reads.out.error_sample
        )

        error_model_ch = learn_errors.out.error_model
    }

    denoise_input_ch = dada2_filter_ccs.out.filtered_fastq
        .combine(error_model_ch)
        .map { sampleID, fastq, error_model ->
            tuple(sampleID, fastq, error_model)
        }

    dada2_denoise_independent(denoise_input_ch)

    dada2_make_seqtab(
        dada2_denoise_independent.out.dada_rds
            .map { sampleID, rds -> rds }
            .collect()
    )
 
    dada2_remove_chimeras(
        dada2_make_seqtab.out.seqtab_rds
    )

    dada2_filter_asvs(
        dada2_remove_chimeras.out.seqtab_nochim_rds,
        params.min_asv_total_freq,
        params.min_asv_sample
    )

    dada2_stats(
        dada2_filter_ccs.out.filter_stats.collect(),
        dada2_denoise_independent.out.denoise_stats.collect(),
        dada2_remove_chimeras.out.seqtab_nochim_rds,
        dada2_filter_asvs.out.seqtab_filtered_rds,
        metadata_ch
    )

    emit:
    asv_fasta     = dada2_filter_asvs.out.asv_fasta
    asv_table_tsv = dada2_filter_asvs.out.asv_table_tsv
}
