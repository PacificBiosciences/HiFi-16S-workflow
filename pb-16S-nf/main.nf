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

  Run the pipeline through usual sample manifest and metadata file used in
  QIIME 2:

  nextflow run main.nf --input samples.tsv --metadata metadata.tsv
  """
}

params.min_len = 1000
params.max_len = 1600
params.front_p = 'AGRGTTYGATYMTGGCTCAG'
params.adapter_p = 'RGYTACCTTGTTACGACTT'
params.pooling_method = 'pseudo'
params.vsearch_db = "~/references/taxonomy_database/qiime2/silva-138-99-seqs.qza"
params.vsearch_tax = "~/references/taxonomy_database/qiime2/silva-138-99-tax.qza"
params.maxreject = 100
params.maxaccept = 5
// TODO rarefaction depth change to max of sample depth dynamically
params.rarefaction_depth = 10000

// Show help message
params.help = false
if (params.help) exit 0, helpMessage()

process sanity_check {
  input:
  path sample_manifest

  output:
  path 'sanity.txt'

  script:
  """
  echo "Running on file $sample_manifest and the file content is:" > sanity.txt
  cat $sample_manifest >> sanity.txt
  """
}

// Import data into QIIME 2
process import_qiime2 {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "import_qiime"

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
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "summary_demux"

  input:
  path samples_qza

  output:
  path "samples.demux.summary.qzv"

  script:
  """
  qiime demux summarize --i-data $samples_qza \
    --o-visualization samples.demux.summary.qzv
  """

}

process dada2_denoise {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "dada2"
  // Use default compute config
  label 'cpu8'

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
    --p-front \'$params.front_p\' \
    --p-adapter \'$params.adapter_p\' \
    --p-n-threads $task.cpus \
    --p-pooling-method \'$params.pooling_method\'
  """
}

// QC summary for dada2
process dada2_qc {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "results"

  input:
  path asv_stats
  path asv_freq
  path metadata

  output:
  path "*"

  script:
  """
  qiime metadata tabulate --m-input-file $asv_stats \
    --o-visualization dada2_stats.qzv

  qiime feature-table summarize --i-table $asv_freq \
    --o-visualization dada2_table.qzv \
    --m-sample-metadata-file $metadata
  """
}

// Rarefaction visualization
process dada2_rarefaction {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "results"

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

process class_tax {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "dada2"
  // Use default compute config
  label 'cpu8'

  input:
  path asv_seq

  output:
  path "taxonomy.vsearch.qza", emit: tax_vsearch
  path "tax_export"

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
  """
}

process barplot {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "results"

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


workflow qiime2 {
  sample_file = channel.fromPath(params.input)
  metadata_file = channel.fromPath(params.metadata)
  sanity_check(sample_file)
  import_qiime2(sample_file)
  demux_summarize(import_qiime2.out)
  dada2_denoise(import_qiime2.out)
  dada2_qc(dada2_denoise.out.asv_stats, dada2_denoise.out.asv_freq, metadata_file)
  dada2_rarefaction(dada2_denoise.out.asv_freq, metadata_file)
  class_tax(dada2_denoise.out.asv_seq)
  barplot(dada2_denoise.out.asv_freq, class_tax.out.tax_vsearch, metadata_file)
}

workflow {
  qiime2()
}
