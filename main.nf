/*
===============================================================================

Nextflow pipeline using QIIME 2 to process CCS data via DADA2 plugin. Takes
in demultiplexed 16S amplicon sequencing FASTQ file.

===============================================================================

Author: Khi Pin, Chua
Updated: 2024-12-17
*/

nextflow.enable.dsl=2
version = "0.9"

// Import modules
include { write_log } from './modules/utils'
include { QC_fastq; QC_fastq_post_trim; collect_QC; collect_QC_skip_cutadapt } from './modules/qc'
include { cutadapt } from './modules/qc'
include { prepare_qiime2_manifest; prepare_qiime2_manifest_skip_cutadapt; merge_manifest } from './modules/qc'
include { import_qiime2; demux_summarize } from './modules/qc'

include { learn_error; dada2_denoise; dada2_denoise_with_error_model } from './modules/dada2'
include { mergeASV; filter_dada2; dada2_qc; dada2_rarefaction } from './modules/dada2'

include { class_tax; dada2_assignTax } from './modules/taxonomy'
include { export_biom; export_biom_skip_nb } from './modules/taxonomy'

include { qiime2_phylogeny_diversity } from './modules/visualization'
include { barplot; barplot_nb; html_rep; html_rep_skip_cutadapt } from './modules/visualization'

include { picrust2; download_db; helpMessage } from './modules/utils'

// Show help message
if (params.help) exit 0, helpMessage()
if (params.version) exit 0, log.info(version)

// Dynamic parameter settings based on input
if (params.input) {
    n_sample = file(params.input).countLines() - 1
    if (n_sample == 1) {
        dynamic_min_asv_totalfreq = 0
        dynamic_min_asv_sample = 0
        println("Only 1 sample. min_asv_sample and min_asv_totalfreq set to 0.")
    } else {
        dynamic_min_asv_totalfreq = params.min_asv_totalfreq
        dynamic_min_asv_sample = params.min_asv_sample
    }
} else {
    println("No input file given to --input!")
    n_sample = 0
    dynamic_min_asv_totalfreq = 0
    dynamic_min_asv_sample = 0
}

// Set primer trimming parameters
if (params.skip_primer_trim) {
    dynamic_front_p = 'none'
    dynamic_adapter_p = 'none'
    trim_cutadapt = "No"
} else {
    dynamic_front_p = params.front_p
    dynamic_adapter_p = params.adapter_p
    trim_cutadapt = "Yes"
}

// Generate log text for parameters
log_text = """
    Parameters set for pb-16S-nf pipeline for PacBio HiFi 16S
    =========================================================
    Number of samples in samples TSV: $n_sample
    Filter input reads above Q: $params.filterQ
    Trim primers with cutadapt: $trim_cutadapt
    Limit to N reads if exceeding N reads (0 = disabled): $params.downsample
    Forward primer: $params.front_p
    Reverse primer: $params.adapter_p
    Minimum amplicon length filtered in DADA2: $params.min_len
    Maximum amplicon length filtered in DADA2: $params.max_len
    Using specified FASTQ for DADA2 error model: $params.learn_error_sample
    DADA2 OMEGA_C parameter (default 1e-40): $params.omegac
    maxEE parameter for DADA2 filterAndTrim: $params.max_ee
    minQ parameter for DADA2 filterAndTrim: $params.minQ
    Pooling method for DADA2 denoise process: $params.pooling_method
    Minimum number of samples required to keep any ASV: $dynamic_min_asv_sample
    Minimum number of reads required to keep any ASV: $dynamic_min_asv_totalfreq 
    Taxonomy sequence database for VSEARCH: $params.vsearch_db
    Taxonomy annotation database for VSEARCH: $params.vsearch_tax
    Skip Naive Bayes classification: $params.skip_nb
    SILVA database for Naive Bayes classifier: $params.silva_db
    GTDB database for Naive Bayes classifier: $params.gtdb_db
    GreenGenes2 database for Naive Bayes classifier: $params.gg2_db
    Database to prioritize for Naive Bayes classifier: $params.db_to_prioritize
    VSEARCH maxreject: $params.maxreject
    VSEARCH maxaccept: $params.maxaccept
    VSEARCH perc-identity: $params.vsearch_identity
    QIIME 2 rarefaction curve sampling depth: $params.rarefaction_depth
    Number of threads specified for cutadapt: $params.cutadapt_cpu
    Number of threads specified for DADA2: $params.dada2_cpu
    Number of threads specified for VSEARCH: $params.vsearch_cpu
    Script location for HTML report generation: $params.rmd_vis_biom_script
    Container enabled via docker/singularity: $params.enable_container
    Version of Nextflow pipeline: $version
"""

workflow pb16S {
    if (params.download_db) {
        download_db()
    } else {
        log.info(log_text)
        write_log(log_text)
        
        // Create input channels
        sample_file = channel
            .fromPath(params.input)
            .splitCsv(header: ['sample', 'fastq'], skip: 1, sep: "\t")
            .map{ row -> tuple(row.sample, file(row.fastq)) }
        
        metadata_file = channel.fromPath(params.metadata)

        // Initial QC
        QC_fastq(sample_file)

        if (params.skip_primer_trim) {
            // Skip primer trimming workflow branch
            collect_QC_skip_cutadapt(
                QC_fastq.out.all_seqkit_stats.collect(),
                QC_fastq.out.all_seqkit_summary.collect()
            )
            
            prepare_qiime2_manifest_skip_cutadapt(QC_fastq.out.filtered_fastq.collect(), metadata_file)
            qiime2_manifest = prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file.flatten()
            merge_manifest(prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file.collect())
            
            
            cutadapt_summary = "none"
            post_trim_readstats = "none"
            collect_QC_readstats = collect_QC_skip_cutadapt.out.all_samples_readstats
            collect_QC_summarised_sample_stats = collect_QC_skip_cutadapt.out.summarised_sample_readstats
            import_qiime2(qiime2_manifest, QC_fastq.out.filtered_fastq_files.collect())
        } else {
            // Normal workflow branch with primer trimming
            cutadapt(
                QC_fastq.out.filtered_fastq,
                dynamic_front_p,
                dynamic_adapter_p
            )
            QC_fastq_post_trim(cutadapt.out.cutadapt_fastq)
            
            collect_QC(
                QC_fastq.out.all_seqkit_stats.collect(),
                QC_fastq.out.all_seqkit_summary.collect(),
                cutadapt.out.summary_tocollect.collect(),
                QC_fastq_post_trim.out.all_seqkit_stats.collect()
            )
            
            prepare_qiime2_manifest(cutadapt.out.cutadapt_fastq.collect(), metadata_file)
            qiime2_manifest = prepare_qiime2_manifest.out.sample_trimmed_file.flatten()
            merge_manifest(prepare_qiime2_manifest.out.sample_trimmed_file.collect())
            
            
            cutadapt_summary = collect_QC.out.cutadapt_summary
            post_trim_readstats = collect_QC.out.all_samples_readstats_post_trim
            collect_QC_readstats = collect_QC.out.all_samples_readstats
            collect_QC_summarised_sample_stats = collect_QC.out.summarised_sample_readstats
            import_qiime2(qiime2_manifest, cutadapt.out.cutadapt_fastq_files.collect())
        }

        // QIIME2 import and demux summary
        demux_summarize(import_qiime2.out)

        // DADA2 processing
        if (params.learn_error_sample) {
            learn_error(params.learn_error_sample, params.learnError_script)
            dada2_denoise_with_error_model(import_qiime2.out, params.dadaCCS_script, params.minQ, learn_error.out.dada2_error_model)
            mergeASV(
                dada2_denoise_with_error_model.out.asv_seq.collect(), 
                dada2_denoise_with_error_model.out.asv_freq.collect(), 
                dada2_denoise_with_error_model.out.asv_stats.collect()
            )
        } else {
            dada2_denoise(import_qiime2.out, params.dadaCCS_script, params.minQ)
            mergeASV(
                dada2_denoise.out.asv_seq.collect(), 
                dada2_denoise.out.asv_freq.collect(), 
                dada2_denoise.out.asv_stats.collect()
            )
        }

        // Filter and QC DADA2 results
        filter_dada2(mergeASV.out.asv_freq, mergeASV.out.asv_seq, dynamic_min_asv_totalfreq, dynamic_min_asv_sample)
        dada2_qc(mergeASV.out.asv_stats, filter_dada2.out.asv_freq, metadata_file)

        // Rarefaction and phylogeny
        if (params.rarefaction_depth > 0) {
            qiime2_phylogeny_diversity(metadata_file, filter_dada2.out.asv_seq,
                filter_dada2.out.asv_freq, params.rarefaction_depth, n_sample)
        } else {
            qiime2_phylogeny_diversity(metadata_file, filter_dada2.out.asv_seq,
                filter_dada2.out.asv_freq, dada2_qc.out.rarefaction_depth, n_sample)
        }
        
        dada2_rarefaction(filter_dada2.out.asv_freq, metadata_file, dada2_qc.out.alpha_depth)

        // Taxonomy classification
        class_tax(filter_dada2.out.asv_seq, filter_dada2.out.asv_freq, 
                 params.vsearch_db, params.vsearch_tax)

        if (params.skip_nb) {
            export_biom_skip_nb(filter_dada2.out.asv_freq, class_tax.out.tax_tsv)
            barplot(filter_dada2.out.asv_freq, class_tax.out.tax_vsearch, 
                   metadata_file, "taxonomy_barplot_vsearch.qzv")
            nb_tax = []
            biom_file = export_biom_skip_nb.out.biom_vsearch
        } else {
            dada2_assignTax(filter_dada2.out.asv_seq_fasta, filter_dada2.out.asv_seq, 
                          filter_dada2.out.asv_freq, params.silva_db, params.gtdb_db, 
                          params.gg2_db, params.dadaAssign_script, params.db_to_prioritize)
            export_biom(filter_dada2.out.asv_freq, dada2_assignTax.out.best_nb_tax, 
                       class_tax.out.tax_tsv)
            barplot(filter_dada2.out.asv_freq, class_tax.out.tax_vsearch, 
                   metadata_file, "taxonomy_barplot_vsearch.qzv")
            barplot_nb(filter_dada2.out.asv_freq, dada2_assignTax.out.best_nb_tax_qza, 
                      metadata_file, "taxonomy_barplot_nb.qzv")
            nb_tax = dada2_assignTax.out.best_nb_tax_tsv
            biom_file = export_biom.out.biom_vsearch
        }

        // Optional PICRUSt2 analysis
        if (params.run_picrust2) {
            picrust2(filter_dada2.out.asv_seq_fasta, biom_file)
        }

        // Generate reports
        if (params.skip_primer_trim) {
            html_rep_skip_cutadapt(
                nb_tax, metadata_file, merge_manifest.out.merged_manifest,
                dada2_qc.out.dada2_qc_tsv, collect_QC_readstats,
                collect_QC_summarised_sample_stats, cutadapt_summary,
                class_tax.out.tax_freq_tab_tsv, qiime2_phylogeny_diversity.out.bray_mat,
                qiime2_phylogeny_diversity.out.unifrac_mat,
                qiime2_phylogeny_diversity.out.wunifrac_mat, params.colorby,
                post_trim_readstats, params.rmd_vis_biom_script, params.rmd_helper
            )
        } else {
            html_rep(
                nb_tax, metadata_file, merge_manifest.out.merged_manifest,
                dada2_qc.out.dada2_qc_tsv, collect_QC_readstats,
                collect_QC_summarised_sample_stats, cutadapt_summary,
                class_tax.out.tax_freq_tab_tsv, qiime2_phylogeny_diversity.out.bray_mat,
                qiime2_phylogeny_diversity.out.unifrac_mat,
                qiime2_phylogeny_diversity.out.wunifrac_mat, params.colorby,
                post_trim_readstats, params.rmd_vis_biom_script, params.rmd_helper
            )
        }
    }
}

workflow {
    pb16S()
}
