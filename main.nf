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

  Sequences can be trimmed first with lima (higher rate compared to using DADA2
  ) using the example command (16S_primers.fasta must be disambiguated first 
  using all possible combinations of degenerate sequences, min-score-lead set 
  to 0 as we don't care if the primers pair are similar, they should be since 
  it's just degenerate sequences!):

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
  --outdir    Output directory name (default "results")
  --vsearch_db	Directory for VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Directory for VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --front_p    Forward 16S primer to trim using DADA2. Set to 'none' if primers already
               trimmed using lima (Default "AGRGTTYGATYMTGGCTCAG")
  --adapter_p    Reverse 16S primer to trim using DADA2. Set to 'none' if primers already
               trimmed using lima (Default "RGYTACCTTGTTACGACTT")
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
params.dada2_cpu = 8
params.vsearch_cpu = 8
params.outdir = "results"
params.max_ee = 2

// Show help message
params.help = false
if (params.help) exit 0, helpMessage()

process sanity_check {
  label 'cpu_def'
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
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "$params.outdir/summary_demux"
  label 'cpu_def'

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
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "$params.outdir/results"
  label 'cpu_def'

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
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "$params.outdir/results"
  cpus params.vsearch_cpu

  input:
  path asv_seq
  path asv_freq

  output:
  path "taxonomy.vsearch.qza", emit: tax_vsearch
  path "tax_export/taxonomy.tsv", emit: tax_tsv
  path "merged_freq_tax.qzv", emit: tax_freq_tab
  path "merged_freq_tax_tsv"

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
  """
}

// Export results into biom for use with phyloseq
process export_biom {
  conda 'qiime2-2022.2-py38-linux-conda.yml'
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path asv_freq
  path tax_tsv

  output:
  path "feature-table-tax.biom"

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
  conda 'qiime2-2022.2-py38-linux-conda.yml'
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

// TODO Export biom
// TODO Export table of ASV to taxonomy and read counts (https://forum.qiime2.org/t/combine-the-taxonomy-table-with-the-asv-count-table/18852/5)
// TODO Visualization of all results in a nice report
// TODO Add QC of input reads. QV, read length etc using seqkits
// TODO Add Krona plot
// TODO Extract Qiime 2 log file under /tmp especially for denoise step
// TODO Classified species percentage
// TODO How to BLAST ASVs for strain level assignment
// TODO Phylogenetic tree. ONT has an interactive one which is nice
// TODO Need to add option to *not* trim primers as MAS-seq may already trimmed away primers
// TODO Look into use Kraken 2

workflow qiime2 {
  sample_file = channel.fromPath(params.input)
  metadata_file = channel.fromPath(params.metadata)
  sanity_check(sample_file)
  import_qiime2(sample_file)
  demux_summarize(import_qiime2.out)
  dada2_denoise(import_qiime2.out)
  dada2_qc(dada2_denoise.out.asv_stats, dada2_denoise.out.asv_freq, metadata_file)
  dada2_rarefaction(dada2_denoise.out.asv_freq, metadata_file)
  class_tax(dada2_denoise.out.asv_seq, dada2_denoise.out.asv_freq)
  export_biom(dada2_denoise.out.asv_freq, class_tax.out.tax_tsv)
  barplot(dada2_denoise.out.asv_freq, class_tax.out.tax_vsearch, metadata_file)
}

workflow {
  qiime2()
}
