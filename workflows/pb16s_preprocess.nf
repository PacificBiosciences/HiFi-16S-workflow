nextflow.enable.dsl = 2

include {
    inspect_metadata
    QC_raw_stats
    filter_fastq
    downsample_fastq
    cutadapt
    summarize_cutadapt
    QC_fastq_post_trim
    collect_QC
    collect_QC_skip_cutadapt
} from '../modules/qc'


workflow PB16S_PREPROCESS {

    take:
    sample_sheet
    metadata

    main:

    inspect_metadata(sample_sheet, metadata)

    sample_ch = sample_sheet
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row['sample-id'] || !row['filepath']) {
                error "Input TSV must contain columns 'sample-id' and 'filepath'"
            }

            tuple(
                row['sample-id'] as String,
                file(row['filepath'] as String)
            )
        }

    QC_raw_stats(sample_ch)
    filter_fastq(sample_ch)

    if (params.skip_primer_trim) {

        downsample_fastq(
            filter_fastq.out.filtered_fastq
        )

        collect_QC_skip_cutadapt(
            QC_raw_stats.out.readstats.collect(),
            QC_raw_stats.out.summarystats.collect()
        )

        reads_for_dada2 = downsample_fastq.out.downsampled_fastq
    }
    else {

        cutadapt(
            filter_fastq.out.filtered_fastq,
            params.forward_p,
            params.reverse_p
        )

        summarize_cutadapt(
            cutadapt.out.cutadapt_report.collect()
        )

        QC_fastq_post_trim(
            cutadapt.out.cutadapt_fastq
        )

        downsample_fastq(
            cutadapt.out.cutadapt_fastq
        )

        collect_QC(
            QC_raw_stats.out.readstats.collect(),
            QC_raw_stats.out.summarystats.collect(),
            summarize_cutadapt.out.summary_tocollect,
            QC_fastq_post_trim.out.all_seqkit_stats.collect()
        )

        reads_for_dada2 = downsample_fastq.out.downsampled_fastq
    }

    emit:
    reads_for_dada2 = reads_for_dada2
}
