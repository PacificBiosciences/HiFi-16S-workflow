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

  nextflow run main.nf --input samples.tsv --metadata metadata.tsv \\
    --dada2_cpu 8 --vsearch_cpu 8

  By default, sequences are first trimmed with cutadapt (higher rate compared to using DADA2
  ) using the example command below. You can skip this by specifying "--skip_primer_trim" 
  if the sequences are already trimmed. The primer sequences used are the F27 and R1492
  primers for full length 16S sequencing.

  Other important options:
  --filterQ    Filter input reads above this Q value (default: 30).
  --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
              this value will be discarded (default: 2)
  --min_len    Minimum length of sequences to keep (default: 1000)
  --max_len    Maximum length of sequences to keep (default: 1600)
  --pooling_method    QIIME 2 pooling method for DADA2 denoise (default: "pseudo"),
                      see QIIME 2 documentation for more details
  --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 5)
  --min_asv_totalfreq    Total frequency of any ASV must be above this threshold
                         across all samples to be retained. Set this to 0 to disable filtering
                         (default 5)
  --min_asv_sample    ASV must exist in at least min_asv_sample to be retained. 
                      Set this to 0 to disable. (default 2)
  --vsearch_identity    Minimum identity to be considered as hit (default 0.97)
  --rarefaction_depth    Rarefaction curve "max-depth" parameter. By default the pipeline
                         automatically select a cut-off above the minimum of the denoised 
                         reads for >90% of the samples. This cut-off is stored in a file called
                         "rarefaction_depth_suggested.txt" file in the results folder
                         (default: null)
  --dada2_cpu    Number of threads for DADA2 denoising (default: 8)
  --vsearch_cpu    Number of threads for VSEARCH taxonomy classification (default: 8)
  --cutadapt_cpu    Number of threads for primer removal using cutadapt (default: 16)
  --outdir    Output directory name (default: "results")
  --vsearch_db	Location of VSEARCH database (e.g. silva-138-99-seqs.qza can be
                downloaded from QIIME database)
  --vsearch_tax    Location of VSEARCH database taxonomy (e.g. silva-138-99-tax.qza can be
                   downloaded from QIIME database)
  --silva_db   Location of Silva 138 database for taxonomy classification 
  --gtdb_db    Location of GTDB r202 for taxonomy classification
  --refseq_db    Location of RefSeq+RDP database for taxonomy classification
  --skip_primer_trim    Skip all primers trimming (switch off cutadapt and DADA2 primers
                        removal) (default: trim with cutadapt)
  --colorby    Columns in metadata TSV file to use for coloring the MDS plot
               in HTML report (default: condition)
  """
}

// Show help message
params.help = false
if (params.help) exit 0, helpMessage()

params.skip_primer_trim = false
params.filterQ = 30
params.min_len = 1000
params.max_len = 1600
params.colorby = "condition"
params.skip_phylotree = false
n_sample=file(params.input).countLines() - 1
if (n_sample == 1) {
  params.min_asv_totalfreq = 0
  params.min_asv_sample = 0
  println("Only 1 sample. min_asv_sample and min_asv_totalfreq set to 0.")
} else {
  params.min_asv_totalfreq = 5
  params.min_asv_sample = 2
}

if (params.skip_primer_trim) {
  params.front_p = 'none'
  params.adapter_p = 'none'
  trim_cutadapt = "No"
} else {
  trim_cutadapt = "Yes"
  params.front_p = 'none'
  params.adapter_p = 'none'
}
// Hidden parameters in case need to trim with DADA2 removePrimers function
// params.front_p = 'AGRGTTYGATYMTGGCTCAG'
// params.adapter_p = 'RGYTACCTTGTTACGACTT'
params.pooling_method = 'pseudo'
params.vsearch_db = "~/references/taxonomy_database/qiime2/silva-138-99-seqs.qza"
params.vsearch_tax = "~/references/taxonomy_database/qiime2/silva-138-99-tax.qza"
params.silva_db = "~/references/taxonomy_database/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
params.gtdb_db = "~/references/taxonomy_database/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz"
params.refseq_db = "~/references/taxonomy_database/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz"
params.maxreject = 100
params.maxaccept = 5
params.rarefaction_depth = null
params.dada2_cpu = 8
params.vsearch_cpu = 8
params.vsearch_identity = 0.97
params.outdir = "results"
params.max_ee = 2
params.rmd_vis_biom_script= "$projectDir/scripts/visualize_biom.Rmd"
// Helper script for biom vis script
params.rmd_helper = "$projectDir/scripts/import_biom.R"
params.cutadapt_cpu = 16
params.primer_fasta = "$projectDir/scripts/16S_primers.fasta"
params.dadaCCS_script = "$projectDir/scripts/run_dada_ccs.R"
params.dadaAssign_script = "$projectDir/scripts/dada2_assign_tax.R"
// kraken2 DB if enabling Kraken2. Deprecated.
// params.kraken2_db = "~/references/taxonomy_database/ncbi_bacteria_16S_k2db/"



log.info """
  Parameters set for pb-16S-nf pipeline for PacBio HiFi 16S
  =========================================================
  Number of samples in samples TSV: $n_sample
  Filter input reads above Q: $params.filterQ
  Trim primers with cutadapt: $trim_cutadapt
  Minimum amplicon length filtered in DADA2: $params.min_len
  Maximum amplicon length filtered in DADA2: $params.max_len
  maxEE parameter for DADA2 filterAndTrim: $params.max_ee
  Pooling method for DADA2 denoise process: $params.pooling_method
  Minimum number of samples required to keep any ASV: $params.min_asv_sample
  Minimum number of reads required to keep any ASV: $params.min_asv_totalfreq 
  Taxonomy sequence database for VSEARCH: $params.vsearch_db
  Taxonomy annotation database for VSEARCH: $params.vsearch_tax
  VSEARCH maxreject: $params.maxreject
  VSEARCH maxaccept: $params.maxaccept
  VSEARCH perc-identity: $params.vsearch_identity
  QIIME 2 rarefaction curve sampling depth: $params.rarefaction_depth
  Number of threads specified for cutadapt: $params.cutadapt_cpu
  Number of threads specified for DADA2: $params.dada2_cpu
  Number of threads specified for VSEARCH: $params.vsearch_cpu
  Script location for HTML report generation: $params.rmd_vis_biom_script
"""

// QC before cutadapt
process QC_fastq {
  conda "$projectDir/pb-16s-pbtools.yml"
  label 'cpu8'
  publishDir "$params.outdir/filtered_input_FASTQ", pattern: '*filterQ*.fastq.gz'

  input:
  tuple val(sampleID), path(sampleFASTQ)

  output:
  path "${sampleID}.seqkit.readstats.tsv", emit: all_seqkit_stats
  path "${sampleID}.seqkit.summarystats.tsv", emit: all_seqkit_summary
  path "${sampleID}_filtered.tsv", emit: filtered_fastq_tsv
  tuple val(sampleID), path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq

  script:
  """
  seqkit fx2tab -j $task.cpus -q --gc -l -H -n $sampleFASTQ |\
    csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.readstats.tsv
  seqkit stats -T -j $task.cpus -a ${sampleFASTQ} |\
    csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.summarystats.tsv
  seqkit seq -j $task.cpus --min-qual $params.filterQ $sampleFASTQ --out-file ${sampleID}.filterQ${params.filterQ}.fastq.gz
  echo -e "${sampleID}\t"\$PWD"/${sampleID}.filterQ${params.filterQ}.fastq.gz" >> ${sampleID}_filtered.tsv
  """
}

// Trim full length 16S primers with cutadapt
process cutadapt {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/trimmed_primers_FASTQ", pattern: '*.fastq.gz'
  publishDir "$params.outdir/cutadapt_summary", pattern: '*.report'
  cpus params.cutadapt_cpu

  input:
  tuple val(sampleID), path(sampleFASTQ)

  output:
  tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq
  path "sample_ind_*.tsv", emit: samples_ind
  path "*.report", emit: cutadapt_summary
  path "cutadapt_summary_${sampleID}.tsv", emit: summary_tocollect

  script:
  """
  cutadapt -a "^AGRGTTYGATYMTGGCTCAG...AAGTCGTAACAAGGTARCY\$" \
    ${sampleFASTQ} \
    -o ${sampleID}.trimmed.fastq.gz \
    -j ${task.cpus} --trimmed-only --revcomp -e 0.1 \
    --json ${sampleID}.cutadapt.report

  echo -e "${sampleID}\t"\$PWD"/${sampleID}.trimmed.fastq.gz" > sample_ind_${sampleID}.tsv
  input_read=`jq -r '.read_counts | .input' ${sampleID}.cutadapt.report`
  demux_read=`jq -r '.read_counts | .output' ${sampleID}.cutadapt.report`
  echo -e "sample\tinput_reads\tdemuxed_reads" > cutadapt_summary_${sampleID}.tsv
  echo -e "${sampleID}\t\$input_read\t\$demux_read" >> cutadapt_summary_${sampleID}.tsv
  """
}

// To process with Kraken2. Deprecated.
// process kraken2 {
//   conda "$projectDir/pb-16s-pbtools.yml"
//   publishDir "$params.outdir/kraken2", pattern: '*kraken2.report'
//   publishDir "$params.outdir/kraken2", pattern: '*kraken2.out'
//   publishDir "$params.outdir/bracken", pattern: '*bracken.report'
//   publishDir "$params.outdir/bracken", pattern: '*bracken.out'
//   label 'cpu8'
// 
//   input:
//   tuple val(sampleID), path(sampleTrimmed_fastq)
// 
//   output:
//   path "${sampleID}.kraken2.report"
//   path "${sampleID}.kraken2.out"
//   path "${sampleID}.bracken.report"
//   path "${sampleID}.bracken.out"
// 
//   script:
//   """
//   kraken2 --db $params.kraken2_db \
//     --report ${sampleID}.kraken2.report \
//     --output ${sampleID}.kraken2.out --use-names \
//     ${sampleTrimmed_fastq} > ${sampleID}.kraken2.log
// 
//   bracken -d $params.kraken2_db -r 1000 -l S -t 0 \
//     -i ${sampleID}.kraken2.report \
//     -o ${sampleID}.bracken.out \
//     -w ${sampleID}.bracken.report > ${sampleID}.bracken.log
//   """
// }

// Collect QC into single files
process collect_QC {
  conda "$projectDir/pb-16s-pbtools.yml"
  publishDir "$params.outdir/results/reads_QC"
  label 'cpu8'

  input:
  path "*"
  path "*"
  path "*"

  output:
  path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
  path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
  path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
  path "seqkit.summarised_stats.group_by_samples.pretty.tsv"
  path "all_samples_cutadapt_stats.tsv", emit: cutadapt_summary

  script:
  """
  csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
  csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv
  csvtk concat -t cutadapt_summary*.tsv > all_samples_cutadapt_stats.tsv
  # Summary read_qual for each sample
  csvtk summary -t -C '%' -g sample -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv
  csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv > seqkit.summarised_stats.group_by_samples.pretty.tsv
  """
}

// Collect QC into single files if skipping cutadapt
process collect_QC_skip_cutadapt {
  conda "$projectDir/pb-16s-pbtools.yml"
  publishDir "$params.outdir/results/reads_QC"
  label 'cpu8'

  input:
  path "*"
  path "*"

  output:
  path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
  path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
  path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
  path "seqkit.summarised_stats.group_by_samples.pretty.tsv"

  script:
  """
  csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
  csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv
  # Summary read_qual for each sample
  csvtk summary -t -C '%' -g sample -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv
  csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv > seqkit.summarised_stats.group_by_samples.pretty.tsv
  """
}

// Put all trimmed samples into a single file for QIIME2 import
process prepare_qiime2_manifest {
  label 'cpu_def'
    publishDir "$params.outdir/results/"

  input: 
  path "sample_ind_*.tsv"
  path "cutadapt_summary_*.tsv"

  output:
  path "samplefile.txt", emit: sample_trimmed_file
  path "samples_demux_rate.tsv"

  """
  echo -e "sample-id\tabsolute-filepath" > samplefile.txt
  cat sample_ind_*.tsv >> samplefile.txt
  cat cutadapt_summary_*.tsv > samples_demux_rate.tsv
  """

}

// Put all trimmed samples into a single file for QIIME2 import if skipping cutadapt
process prepare_qiime2_manifest_skip_cutadapt {
  label 'cpu_def'
    publishDir "$params.outdir/results/"

  input: 
  path "*filtered.tsv"

  output:
  path "samplefile.txt", emit: sample_trimmed_file

  """
  echo -e "sample-id\tabsolute-filepath" > samplefile.txt
  cat *filtered.tsv >> samplefile.txt
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
  path "dada2-ccs_rep_filtered.qza", emit: asv_seq
  path "dada2-ccs_table_filtered.qza", emit: asv_freq
  path "dada2-ccs_stats.qza", emit:asv_stats
  path "dada2_ASV.fasta", emit:asv_seq_fasta
  path "seqtab_nochim.rds", emit: dada2_rds

  script:
  """
  # Use custom script that can skip primer trimming
  cp $params.dadaCCS_script run_dada_ccs.R
  chmod +x run_dada_ccs.R
  export PATH="./:\$PATH"
  which run_dada_ccs.R
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

  if [ $params.min_asv_sample -gt 0 ] && [ $params.min_asv_totalfreq -gt 0 ]
  then
    # Filter ASVs and sequences
    qiime feature-table filter-features \
      --i-table dada2-ccs_table.qza \
      --p-min-frequency $params.min_asv_totalfreq \
      --p-min-samples $params.min_asv_sample \
      --p-filter-empty-samples \
      --o-filtered-table dada2-ccs_table_filtered.qza
  elif [ $params.min_asv_sample -gt 0 ] && [ $params.min_asv_totalfreq -eq 0 ]
  then
    qiime feature-table filter-features \
      --i-table dada2-ccs_table.qza \
      --p-min-samples $params.min_asv_sample \
      --p-filter-empty-samples \
      --o-filtered-table dada2-ccs_table_filtered.qza
  elif [ $params.min_asv_sample -eq 0 ] && [ $params.min_asv_totalfreq -gt 0 ]
  then
    qiime feature-table filter-features \
      --i-table dada2-ccs_table.qza \
      --p-min-frequency $params.min_asv_totalfreq \
      --p-filter-empty-samples \
      --o-filtered-table dada2-ccs_table_filtered.qza
  else
    mv dada2-ccs_table.qza dada2-ccs_table_filtered.qza
  fi
  qiime feature-table filter-seqs \
    --i-data dada2-ccs_rep.qza \
    --i-table dada2-ccs_table_filtered.qza \
    --o-filtered-data dada2-ccs_rep_filtered.qza

  # Export FASTA file for ASVs
  qiime tools export --input-path dada2-ccs_rep_filtered.qza \
    --output-path .
  mv dna-sequences.fasta dada2_ASV.fasta
  """
}


// Assign taxonomies to SILVA, GTDB and RefSeq using DADA2
// assignTaxonomy function based on Naive Bayes classifier
process dada2_assignTax {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  cpus params.vsearch_cpu

  input:
  path asv_seq_fasta
  path asv_seq
  path asv_freq

  output:
  path "best_tax.qza", emit:best_nb_tax_qza
  path "best_taxonomy.tsv", emit: best_nb_tax
  path "best_tax_merged_freq_tax.tsv", emit: best_nb_tax_tsv

  script:
  """
  Rscript --vanilla $params.dadaAssign_script $asv_seq_fasta $task.cpus $params.silva_db $params.gtdb_db $params.refseq_db 80

  qiime feature-table transpose --i-table $asv_freq \
    --o-transposed-feature-table transposed-asv.qza

  qiime tools import --type "FeatureData[Taxonomy]" \
    --input-format "TSVTaxonomyFormat" \
    --input-path best_taxonomy.tsv --output-path best_tax.qza

  qiime metadata tabulate --m-input-file $asv_seq \
    --m-input-file best_tax.qza \
    --m-input-file transposed-asv.qza \
    --o-visualization merged_freq_tax.qzv

  qiime tools export --input-path merged_freq_tax.qzv \
    --output-path merged_freq_tax_tsv

  mv merged_freq_tax_tsv/metadata.tsv best_tax_merged_freq_tax.tsv
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
  env(int_rarefaction_d), emit: rarefaction_depth
  path "rarefaction_depth_suggested.txt"

  script:
  """
  qiime metadata tabulate --m-input-file $asv_stats \
    --o-visualization dada2_stats.qzv

  qiime feature-table summarize --i-table $asv_freq \
    --o-visualization dada2_table.qzv \
    --m-sample-metadata-file $metadata

  qiime tools export --input-path dada2_table.qzv \
    --output-path dada2_table_summary

  qiime tools export --input-path $asv_stats \
    --output-path ./

  # Get number of reads for ASV covering 90% of samples
  number=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | wc -l`
  ninety=0.9
  # Handle single sample and 2 samples
  if [ \${number} == 1 ];
  then
    result=1
  elif [ \${number} == 2 ];
  then
    result=2
  else
    result=`echo "(\$number * \$ninety)" | bc -l` 
  fi
  rarefaction_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
    head -n \${result%%.*} | tail -n1 | cut -f2 -d,`
  int_rarefaction_d=\${rarefaction_d%%.*}
  echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt

  qiime tools export --input-path dada2_stats.qzv \
    --output-path ./dada2_stats
  mv dada2_stats/metadata.tsv dada2_qc.tsv
  """
}

process qiime2_phylogeny_diversity {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results/phylogeny_diversity"
  label 'cpu8' 

  input:
  path metadata
  path asv_seq
  path asv_freq
  val(rarefaction_depth)

  output:
  path "*.qza"
  path "core-metrics-diversity"
  path "phylotree_mafft_rooted.nwk", emit:qiime2_tree
  path "core-metrics-diversity/bray_curtis_distance_matrix.tsv", emit:bray_mat
  path "core-metrics-diversity/weighted_unifrac_distance_matrix.tsv", emit:wunifrac_mat
  path "core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv", emit:unifrac_mat

  script:
  if (n_sample > 1)
  """
  qiime phylogeny align-to-tree-mafft-fasttree \
    --p-n-threads $task.cpus \
    --i-sequences $asv_seq \
    --o-alignment mafft_alignment.qza \
    --o-masked-alignment mafft_alignment_masked.qza \
    --o-tree phylotree_mafft_unrooted.qza \
    --o-rooted-tree phylotree_mafft_rooted.qza

  qiime tools export --input-path phylotree_mafft_rooted.qza \
    --output-path ./
  mv tree.nwk phylotree_mafft_rooted.nwk

  qiime diversity core-metrics-phylogenetic \
    --p-n-jobs-or-threads $task.cpus \
    --i-phylogeny phylotree_mafft_rooted.qza \
    --i-table $asv_freq \
    --m-metadata-file $metadata \
    --p-sampling-depth $rarefaction_depth \
    --output-dir ./core-metrics-diversity

  # Export various matrix for plotting later
  qiime tools export --input-path ./core-metrics-diversity/bray_curtis_distance_matrix.qza \
    --output-path ./core-metrics-diversity
  mv ./core-metrics-diversity/distance-matrix.tsv \
    ./core-metrics-diversity/bray_curtis_distance_matrix.tsv
  qiime tools export --input-path ./core-metrics-diversity/weighted_unifrac_distance_matrix.qza \
    --output-path ./core-metrics-diversity
  mv ./core-metrics-diversity/distance-matrix.tsv \
    ./core-metrics-diversity/weighted_unifrac_distance_matrix.tsv
  qiime tools export --input-path ./core-metrics-diversity/unweighted_unifrac_distance_matrix.qza \
    --output-path ./core-metrics-diversity
  mv ./core-metrics-diversity/distance-matrix.tsv \
    ./core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv
  """
  else
  """
  qiime phylogeny align-to-tree-mafft-fasttree \
    --p-n-threads $task.cpus \
    --i-sequences $asv_seq \
    --o-alignment mafft_alignment.qza \
    --o-masked-alignment mafft_alignment_masked.qza \
    --o-tree phylotree_mafft_unrooted.qza \
    --o-rooted-tree phylotree_mafft_rooted.qza

  qiime tools export --input-path phylotree_mafft_rooted.qza \
    --output-path ./
  mv tree.nwk phylotree_mafft_rooted.nwk
  mkdir ./core-metrics-diversity
  touch ./core-metrics-diversity/bray_curtis_distance_matrix.tsv \
    ./core-metrics-diversity/weighted_unifrac_distance_matrix.tsv \
    ./core-metrics-diversity/unweighted_unifrac_distance_matrix.tsv
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
  val(rarefaction_depth)

  output:
  path "*"

  script:
  if( !params.rarefaction_depth )
    """
    qiime diversity alpha-rarefaction --i-table $asv_freq \
      --m-metadata-file $metadata \
      --o-visualization alpha-rarefaction-curves.qzv \
      --p-min-depth 10 --p-max-depth $rarefaction_depth
    """
  else
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
  path "vsearch_merged_freq_tax.tsv", emit: tax_freq_tab_tsv

  script:
  """
  qiime feature-classifier classify-consensus-vsearch --i-query $asv_seq \
    --o-classification taxonomy.vsearch.qza \
    --i-reference-reads $params.vsearch_db \
    --i-reference-taxonomy $params.vsearch_tax \
    --p-threads $task.cpus \
    --p-maxrejects $params.maxreject \
    --p-maxaccepts $params.maxaccept \
    --p-perc-identity $params.vsearch_identity

  qiime tools export --input-path taxonomy.vsearch.qza --output-path tax_export

  qiime feature-table transpose --i-table $asv_freq \
    --o-transposed-feature-table transposed-asv.qza

  qiime metadata tabulate --m-input-file $asv_seq \
    --m-input-file taxonomy.vsearch.qza \
    --m-input-file transposed-asv.qza \
    --o-visualization merged_freq_tax.qzv

  qiime tools export --input-path merged_freq_tax.qzv \
    --output-path merged_freq_tax_tsv

  mv merged_freq_tax_tsv/metadata.tsv vsearch_merged_freq_tax.tsv
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
  path reads_qc
  path summarised_reads_qc
  path cutadapt_summary_qc
  path vsearch_tax_tsv
  path bray_mat
  path unifrac_mat
  path wunifrac_mat
  val(colorby)

  output:
  path "visualize_biom.html", emit: html_report

  script:
  """
  cp $params.rmd_vis_biom_script visualize_biom.Rmd
  cp $params.rmd_helper import_biom.R
  Rscript -e 'rmarkdown::render("visualize_biom.Rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat"), output_dir="./")'
  """
}

// HTML report
process html_rep_skip_cutadapt {
  conda "$projectDir/pb-16s-vis-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'

  input:
  path tax_freq_tab_tsv
  path metadata
  path sample_manifest
  path dada2_qc
  path reads_qc
  path summarised_reads_qc
  val(cutadapt_summary_qc)
  path vsearch_tax_tsv
  path bray_mat
  path unifrac_mat
  path wunifrac_mat
  val(colorby)

  output:
  path "visualize_biom.html", emit: html_report

  script:
  """
  cp $params.rmd_vis_biom_script visualize_biom.Rmd
  cp $params.rmd_helper import_biom.R
  Rscript -e 'rmarkdown::render("visualize_biom.Rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat"), output_dir="./")'
  """
}

// Krona plot
process krona_plot {
  conda "$projectDir/qiime2-2022.2-py38-linux-conda.yml"
  publishDir "$params.outdir/results"
  label 'cpu_def'
  // Ignore if this fail
  errorStrategy = 'ignore'

  input:
  path asv_freq
  path taxonomy

  output:
  path "krona.qzv"
  path "krona_html"

  script:
  """
  pip install git+https://github.com/kaanb93/q2-krona.git
  qiime krona collapse-and-plot --i-table $asv_freq --i-taxonomy $taxonomy --o-krona-plot krona.qzv
  qiime tools export --input-path krona.qzv --output-path krona_html
  """

}

workflow pb16S {
  sample_file = channel.fromPath(params.input)
    .splitCsv(header: ['sample', 'fastq'], skip: 1, sep: "\t")
    .map{ row -> tuple(row.sample, file(row.fastq)) }
  metadata_file = channel.fromPath(params.metadata)
  QC_fastq(sample_file)
  if (params.skip_primer_trim){
    // kraken2(QC_fastq.out.filtered_fastq)
    collect_QC_skip_cutadapt(QC_fastq.out.all_seqkit_stats.collect(),
        QC_fastq.out.all_seqkit_summary.collect())
    prepare_qiime2_manifest_skip_cutadapt(QC_fastq.out.filtered_fastq_tsv.collect())
    cutadapt_summary = "none"
    qiime2_manifest = prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file
    collect_QC_readstats = collect_QC_skip_cutadapt.out.all_samples_readstats
    collect_QC_summarised_sample_stats = collect_QC_skip_cutadapt.out.summarised_sample_readstats
  } else {
    cutadapt(QC_fastq.out.filtered_fastq)
    // kraken2(cutadapt.out.cutadapt_fastq)
    collect_QC(QC_fastq.out.all_seqkit_stats.collect(),
        QC_fastq.out.all_seqkit_summary.collect(),
        cutadapt.out.summary_tocollect.collect())
    prepare_qiime2_manifest(cutadapt.out.samples_ind.collect(),
        cutadapt.out.summary_tocollect.collect())
    cutadapt_summary = collect_QC.out.cutadapt_summary
    qiime2_manifest = prepare_qiime2_manifest.out.sample_trimmed_file
    collect_QC_readstats = collect_QC.out.all_samples_readstats
    collect_QC_summarised_sample_stats = collect_QC.out.summarised_sample_readstats
  }
  import_qiime2(qiime2_manifest)
  demux_summarize(import_qiime2.out)
  dada2_denoise(import_qiime2.out)
  dada2_qc(dada2_denoise.out.asv_stats, dada2_denoise.out.asv_freq, metadata_file)
  qiime2_phylogeny_diversity(metadata_file, dada2_denoise.out.asv_seq,
      dada2_denoise.out.asv_freq, dada2_qc.out.rarefaction_depth)
  dada2_rarefaction(dada2_denoise.out.asv_freq, metadata_file, dada2_qc.out.rarefaction_depth)
  class_tax(dada2_denoise.out.asv_seq, dada2_denoise.out.asv_freq)
  dada2_assignTax(dada2_denoise.out.asv_seq_fasta, dada2_denoise.out.asv_seq, dada2_denoise.out.asv_freq)
  export_biom(dada2_denoise.out.asv_freq, dada2_assignTax.out.best_nb_tax)
  barplot(dada2_denoise.out.asv_freq, dada2_assignTax.out.best_nb_tax_qza, metadata_file)
  if (params.skip_primer_trim){
    html_rep_skip_cutadapt(dada2_assignTax.out.best_nb_tax_tsv, metadata_file, qiime2_manifest,
        dada2_qc.out.dada2_qc_tsv, 
        collect_QC_readstats, collect_QC_summarised_sample_stats,
        cutadapt_summary, class_tax.out.tax_freq_tab_tsv, qiime2_phylogeny_diversity.out.bray_mat,
        qiime2_phylogeny_diversity.out.unifrac_mat, qiime2_phylogeny_diversity.out.wunifrac_mat,
        params.colorby)
  } else {
    html_rep(dada2_assignTax.out.best_nb_tax_tsv, metadata_file, qiime2_manifest,
        dada2_qc.out.dada2_qc_tsv, 
        collect_QC_readstats, collect_QC_summarised_sample_stats,
        cutadapt_summary, class_tax.out.tax_freq_tab_tsv, qiime2_phylogeny_diversity.out.bray_mat,
        qiime2_phylogeny_diversity.out.unifrac_mat, qiime2_phylogeny_diversity.out.wunifrac_mat,
        params.colorby)
  }
  krona_plot(dada2_denoise.out.asv_freq, dada2_assignTax.out.best_nb_tax_qza)
}

workflow {
  pb16S()
}
