nextflow.enable.dsl = 2

include {
    inspect_metadata
    QC_fastq
    cutadapt
    summarize_cutadapt
    QC_fastq_post_trim
    collect_QC
    collect_QC_skip_cutadapt
} from '../modules/qc'

include { DADA2_WORKFLOW } from './dada2'
include { TAXONOMY_WORKFLOW } from './taxonomy'

workflow PB16S_PREPROCESS {

    take:
    db_manifest

    main:

    sample_sheet_ch = Channel.fromPath(params.input)
    metadata_ch     = Channel.fromPath(params.metadata)

    inspect_metadata(sample_sheet_ch, metadata_ch)

    sample_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row['sample-id'] || !row['filepath']) {
                error "Input TSV must contain columns 'sample-id' and 'filepath'"
            }
            tuple(row['sample-id'] as String, file(row['filepath'] as String))
        }

    QC_fastq(sample_ch)

    if (params.skip_primer_trim) {
        collect_QC_skip_cutadapt(
            QC_fastq.out.all_seqkit_stats.collect(),
            QC_fastq.out.all_seqkit_summary.collect()
        )

        reads_for_dada2 = QC_fastq.out.filtered_fastq
    }
    else {
        cutadapt(
            QC_fastq.out.filtered_fastq,
            params.front_p,
            params.adapter_p
        )

        summarize_cutadapt(cutadapt.out.cutadapt_report)
        QC_fastq_post_trim(cutadapt.out.cutadapt_fastq)

        collect_QC(
            QC_fastq.out.all_seqkit_stats.collect(),
            QC_fastq.out.all_seqkit_summary.collect(),
            summarize_cutadapt.out.summary_tocollect.collect(),
            QC_fastq_post_trim.out.all_seqkit_stats.collect()
        )

        reads_for_dada2 = cutadapt.out.cutadapt_fastq
    }

    DADA2_WORKFLOW(reads_for_dada2, metadata_ch)

    TAXONOMY_WORKFLOW(
        db_manifest,
        DADA2_WORKFLOW.out.asv_fasta,
        DADA2_WORKFLOW.out.asv_table_tsv
    )
}
