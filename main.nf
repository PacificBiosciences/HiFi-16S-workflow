/*
===============================================================================

Nextflow pipeline using QIIME 2 to process CCS data via DADA2 plugin. Takes
in demultiplexed 16S amplicon sequencing FASTQ file.

===============================================================================

Author: Khi Pin, Chua
Updated: 2022-5-5
*/

nextflow.enable.dsl=2

def helpMessage() {
  log.info"""
  Usage:
  This pipeline takes in the standard sample manifest and metadata file used in
  QIIME 2 and produces QC summary, taxonomy classification results and visualization.

  For samples TSV, two columns named "sample-id" and "absolute-filepath" are
  required. For metadata TSV file, at least two columns named "sample_name" and
  "condition" to separate samples into different groups.

  // TODO If no metadata TSV produced, generate a fake one

  nextflow run main.nf --input samples.tsv --metadata metadata.tsv \\
    --dada2_cpu 8 --vsearch_cpu 8

  By default, sequences are first trimmed with lima (higher rate compared to using DADA2
  ) using the example command below. 16S_primers.fasta file contains disambiguated 16S primers
  using all possible combinations of degenerate sequences. "min-score-lead" is set
  to 0 as we don't care if the primers pair are similar, they should be since
  it's just degenerate sequences!

  lima --hifi-preset ASYMMETRIC \\
    demultiplex.16S_For_bc1005--16S_Rev_bc1057.hifi_reads.fastq.gz \\
    16S_primers.fasta \\
    bc1005-bc1057.16s.lima.same.fastq.gz \\
    --log-level INFO \\
    --min-score-lead 0

  Other important options:
  --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
              this value will be discarded (defaul 2)
  --min_len    Minimum length of sequences to keep (default 1000)
  --max_len    Maximum length of sequences to keep (default 1600)
  --pooling_method    QIIME 2 pooling method for DADA2 denoise (default "pseudo"),
                      see QIIME 2 documentation for more details
  --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default 100)
  --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default 5)
  --dada2_cpu    Number of threads for DADA2 denoising (default 8)
  --vsearch_cpu    Number of threads for VSEARCH taxonomy classification (default 8)
  --lima_cpu    Number of threads for primer removal using lima (default 16)
  --outdir    Output directory name (default "results")
  --vsearch_db	Directory for VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Directory for VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --front_p    Forward 16S primer to trim using DADA2. Set to 'none' if primers already
               trimmed using lima (Default "none")
  --adapter_p    Reverse 16S primer to trim using DADA2. Set to 'none' if primers already
               trimmed using lima (Default "none")
  """
}

params.min_len = 1000
params.max_len = 1600
// params.front_p = 'AGRGTTYGATYMTGGCTCAG'
// params.adapter_p = 'RGYTACCTTGTTACGACTT'
params.front_p = 'none'
params.adapter_p = 'none'
params.pooling_method = 'pseudo'
params.vsearch_db = "~/references/taxonomy_database/qiime2/silva-138-99-seqs.qza"
params.vsearch_tax = "~/references/taxonomy_database/qiime2/silva-138-99-tax.qza"
params.maxreject = 100
params.maxaccept = 5
// TODO rarefaction depth change to max of sample depth dynamically. Basically using 
// csvtk to sort and get the count that covers at least 90% of samples
params.rarefaction_depth = 10000
params.dada2_cpu = 8
params.vsearch_cpu = 8
params.outdir = "results"
params.max_ee = 2
params.rmd_vis_biom_script= "$projectDir/scripts/visualize_biom.Rmd"
// Helper script for biom vis script
params.rmd_helper = "$projectDir/scripts/import_biom.R"
params.lima_cpu = 16
params.primer_fasta = "$projectDir/scripts/16S_primers.fasta"

log.info """
  Running pb-16S-nf pipeline for PacBio HiFi 16S
  ==============================================
  Minimum amplicon length filtered in DADA2: $params.min_len
  Maximum amplicon length filtered in DADA2: $params.max_len
  Forward 16S primer trimmed by DADA2 if given: $params.front_p
  Reverse 16S primer trimmed by DADA2 if given: $params.adapter_p
  maxEE parameter for DADA2 filterAndTrim: $params.max_ee
  Pooling method for DADA2 denoise process: $params.pooling_method
  Taxonomy sequence database for VSEARCH: $params.vsearch_db
  Taxonomy annotation database for VSEARCH: $params.vsearch_tax
  VSEARCH maxreject: $params.maxreject
  VSEARCH maxaccept: $params.maxaccept
  QIIME 2 rarefaction curve sampling depth: $params.rarefaction_depth
  Number of threads specified for lima: $params.lima_cpu
  Number of threads specified for DADA2: $params.dada2_cpu
  Number of threads specified for VSEARCH: $params.vsearch_cpu
  Script location for HTML report generation: $params.rmd_vis_biom_script
"""

// Show help message
params.help = false
if (params.help) exit 0, helpMessage()

// Trim full length 16S primers with lima
process lima {
  conda "$projectDir/pb-16s-pbtools.yml"
  publishDir "$params.outdir/trimmed_primers_FASTQ", pattern: '*.fastq.gz'
  publishDir "$params.outdir/lima_summary", pattern: '*.summary'
  cpus params.lima_cpu

  input:
  tuple val(sampleID), path(sampleFASTQ)

  output:
  tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz")
  path "sample_ind_*.tsv", emit: samples_ind
  path "*.summary", emit: lima_summary
  path "lima_summary_${sampleID}.tsv", emit: summary_tocollect

  script:
  """
  lima --hifi-preset ASYMMETRIC --min-score-lead 0 $sampleFASTQ $params.primer_fasta ${sampleID}.trimmed.fastq.gz \
    --log-level INFO --log-file ${sampleID}.lima.log
  echo -e "${sampleID}\t"\$PWD"/${sampleID}.trimmed.fastq.gz" > sample_ind_${sampleID}.tsv
  demux_rate=`grep "ZMWs above all thresholds" *.summary | sed 's/.*(\\([0-9].*%\\))/\\1/g'`
  echo -e "${sampleID}\t\$demux_rate" > lima_summary_${sampleID}.tsv
  """
}

process prepare_qiime2_manifest {
  label 'cpu_def'
    publishDir "$params.outdir/results/"

  input: 
  path "sample_ind_*.tsv"
  path "lima_summary_*.tsv"

  output:
  path "samplefile.txt", emit: sample_trimmed_file
  path "samples_demux_rate.tsv"

  """
  echo -e "sample-id\tabsolute-filepath" > samplefile.txt
  cat sample_ind_*.tsv >> samplefile.txt
  cat lima_summary_*.tsv > samples_demux_rate.tsv
  """

}

// Import data into QIIME 2
process import_qiime2 {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/import_qiime"
  label 'cpu_def'

  input:
  path sample_manifest

  output:
  path 'samples.qza'

  script:
  """
  qiime tools import --type 'SampleData[SequencesWithQuality]' \
    --input-path $sample_manifest \
    --output-path samples.qza \
    --input-format SingleEndFastqManifestPhred33V2
  """
}

process demux_summarize {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/summary_demux"
  label 'cpu_def'

  input:
  path samples_qza

  output:
  path "samples.demux.summary.qzv"
  path "per-sample-fastq-counts.tsv"

  script:
  """
  qiime demux summarize --i-data $samples_qza \
    --o-visualization samples.demux.summary.qzv

  qiime tools export --input-path samples.demux.summary.qzv \
    --output-path ./
  """

}

process dada2_denoise {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/dada2"
  cpus params.dada2_cpu

  input:
  path samples_qza

  output:
  path "dada2-ccs_rep.qza", emit: asv_seq
  path "dada2-ccs_table.qza", emit: asv_freq
  path "dada2-ccs_stats.qza", emit:asv_stats

  script:
  """
  qiime dada2 denoise-ccs --i-demultiplexed-seqs $samples_qza \
    --o-table dada2-ccs_table.qza \
    --o-representative-sequences dada2-ccs_rep.qza \
    --o-denoising-stats dada2-ccs_stats.qza \
    --p-min-len $params.min_len --p-max-len $params.max_len \
    --p-max-ee $params.max_ee \
    --p-front \'$params.front_p\' \
    --p-adapter \'$params.adapter_p\' \
    --p-n-threads $task.cpus \
    --p-pooling-method \'$params.pooling_method\'
  """
}

// QC summary for dada2
process dada2_qc {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path asv_stats
  path asv_freq
  path metadata

  output:
  path "dada2_stats.qzv", emit: dada2_stats
  path "dada2_table.qzv", emit: dada2_table
  path "stats.tsv", emit: dada2_stats_tsv
  path "dada2_qc.tsv", emit: dada2_qc_tsv

  script:
  """
  qiime metadata tabulate --m-input-file $asv_stats \
    --o-visualization dada2_stats.qzv

  qiime feature-table summarize --i-table $asv_freq \
    --o-visualization dada2_table.qzv \
    --m-sample-metadata-file $metadata

  qiime tools export --input-path $asv_stats \
    --output-path ./

  qiime tools export --input-path dada2_stats.qzv \
    --output-path ./dada2_stats
  mv dada2_stats/metadata.tsv dada2_qc.tsv
  """
}

// Rarefaction visualization
process dada2_rarefaction {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path asv_freq
  path metadata

  output:
  path "*"

  script:
  """
  qiime diversity alpha-rarefaction --i-table $asv_freq \
    --m-metadata-file $metadata \
    --o-visualization alpha-rarefaction-curves.qzv \
    --p-min-depth 10 --p-max-depth $params.rarefaction_depth
  """
}

// Classify taxonomy and export table
process class_tax {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  cpus params.vsearch_cpu

  input:
  path asv_seq
  path asv_freq

  output:
  path "taxonomy.vsearch.qza", emit: tax_vsearch
  path "tax_export/taxonomy.tsv", emit: tax_tsv
  path "merged_freq_tax.qzv", emit: tax_freq_tab
  path "merged_freq_tax.tsv", emit: tax_freq_tab_tsv

  script:
  """
  qiime feature-classifier classify-consensus-vsearch --i-query $asv_seq \
    --o-classification taxonomy.vsearch.qza \
    --i-reference-reads $params.vsearch_db \
    --i-reference-taxonomy $params.vsearch_tax \
    --p-threads $task.cpus \
    --p-maxrejects $params.maxreject \
    --p-maxaccepts $params.maxaccept

  qiime tools export --input-path taxonomy.vsearch.qza --output-path tax_export

  qiime feature-table transpose --i-table $asv_freq \
    --o-transposed-feature-table transposed-asv.qza

  qiime metadata tabulate --m-input-file $asv_seq \
    --m-input-file taxonomy.vsearch.qza \
    --m-input-file transposed-asv.qza \
    --o-visualization merged_freq_tax.qzv

  qiime tools export --input-path merged_freq_tax.qzv \
    --output-path merged_freq_tax_tsv

  mv merged_freq_tax_tsv/metadata.tsv merged_freq_tax.tsv
  """
}

// Export results into biom for use with phyloseq
process export_biom {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path asv_freq
  path tax_tsv

  output:
  path "feature-table-tax.biom", emit:biom

  script:
  """
  qiime tools export --input-path $asv_freq --output-path asv_freq/

  sed 's/Feature ID/#OTUID/' $tax_tsv | sed 's/Taxon/taxonomy/' | \
    sed 's/Consensus/confidence/' > biom-taxonomy.tsv

  biom add-metadata -i asv_freq/feature-table.biom \
    -o feature-table-tax.biom \
    --observation-metadata-fp biom-taxonomy.tsv \
    --sc-separated taxonomy
  """
}

process barplot {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path asv_tab
  path tax
  path metadata

  output:
  path 'taxa_barplot.qzv'

  script:
  """
  qiime taxa barplot --i-table $asv_tab --i-taxonomy $tax \
    --m-metadata-file $metadata \
    --o-visualization taxa_barplot.qzv
  """
}

// HTML report
process html_rep {
  conda "$projectDir/pb-16s-vis-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path tax_freq_tab_tsv
  path metadata
  path sample_manifest
  path dada2_qc

  output:
  path "visualize_biom.html", emit: html_report

  script:
  """
  cp $params.rmd_vis_biom_script visualize_biom.Rmd
  cp $params.rmd_helper import_biom.R
  Rscript -e 'rmarkdown::render("visualize_biom.Rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc"), output_dir="./")'
  """

}

// TODO Visualization of all results in a nice report
// TODO Add QC of input reads. QV, read length etc using seqkits
// TODO Add Krona plot
// TODO Extract Qiime 2 log file under /tmp especially for denoise step
// TODO How to BLAST ASVs for strain level assignment
// TODO Phylogenetic tree. ONT has an interactive one which is nice
// TODO Dockerize. Remember to replace run_dada_ccs.R so it works with trimmed sequences
// TODO Look into use Kraken 2

workflow qiime2 {
  sample_file = channel.fromPath(params.input)
    .splitCsv(header: ['sample', 'fastq'], skip: 1, sep: "\t")
    .map{ row -> tuple(row.sample, file(row.fastq)) }
  metadata_file = channel.fromPath(params.metadata)
  lima(sample_file)
  prepare_qiime2_manifest(lima.out.samples_ind.collect(), lima.out.summary_tocollect.collect())
  import_qiime2(prepare_qiime2_manifest.out.sample_trimmed_file)
  demux_summarize(import_qiime2.out)
  dada2_denoise(import_qiime2.out)
  dada2_qc(dada2_denoise.out.asv_stats, dada2_denoise.out.asv_freq, metadata_file)
  dada2_rarefaction(dada2_denoise.out.asv_freq, metadata_file)
  class_tax(dada2_denoise.out.asv_seq, dada2_denoise.out.asv_freq)
  export_biom(dada2_denoise.out.asv_freq, class_tax.out.tax_tsv)
  barplot(dada2_denoise.out.asv_freq, class_tax.out.tax_vsearch, metadata_file)
  html_rep(class_tax.out.tax_freq_tab_tsv, metadata_file, prepare_qiime2_manifest.out.sample_trimmed_file, dada2_qc.out.dada2_qc_tsv)
}

workflow {
  qiime2()
}
