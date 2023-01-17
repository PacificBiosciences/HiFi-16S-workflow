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

  By default, sequences are first trimmed with cutadapt. If adapters are already trimmed, you can skip 
  cutadapt by specifying "--skip_primer_trim".

  Other important options:
  --front_p    Forward primer sequence. Default to F27. (default: AGRGTTYGATYMTGGCTCAG)
  --adapter_p    Reverse primer sequence. Default to R1492. (default: AAGTCGTAACAAGGTARCY)
  --filterQ    Filter input reads above this Q value (default: 20).
  --max_ee    DADA2 max_EE parameter. Reads with number of expected errors higher than
              this value will be discarded (default: 2)
  --minQ    DADA2 minQ parameter. Reads with any base lower than this score 
            will be removed (default: 0)
  --min_len    Minimum length of sequences to keep (default: 1000)
  --max_len    Maximum length of sequences to keep (default: 1600)
  --pooling_method    QIIME 2 pooling method for DADA2 denoise see QIIME 2 
                      documentation for more details (default: "pseudo", alternative: "independent") 
  --maxreject    max-reject parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --maxaccept    max-accept parameter for VSEARCH taxonomy classification method in QIIME 2
                 (default: 100)
  --min_asv_totalfreq    Total frequency of any ASV must be above this threshold
                         across all samples to be retained. Set this to 0 to disable filtering
                         (default 5)
  --min_asv_sample    ASV must exist in at least min_asv_sample to be retained. 
                      Set this to 0 to disable. (default 1)
  --vsearch_identity    Minimum identity to be considered as hit (default 0.97)
  --rarefaction_depth    Rarefaction curve "max-depth" parameter. By default the pipeline
                         automatically select a cut-off above the minimum of the denoised 
                         reads for >80% of the samples. This cut-off is stored in a file called
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
  --skip_nb    Skip Naive-Bayes classification (only uses VSEARCH) (default: false)
  --colorby    Columns in metadata TSV file to use for coloring the MDS plot
               in HTML report (default: condition)
  --run_picrust2    Run PICRUSt2 pipeline. Note that pathway inference with 16S using PICRUSt2
                    has not been tested systematically (default: false)
  --download_db    Download databases needed for taxonomy classification only. Will not
                   run the pipeline. Databases will be downloaded to a folder "databases"
                   in the Nextflow pipeline directory.
  --publish_dir_mode    Outputs mode based on Nextflow "publishDir" directive. Specify "copy"
                        if requires hard copies. (default: symlink)
  --version    Output version
  """
}

// Show help message
params.help = false
if (params.help) exit 0, helpMessage()
params.version = false
version = "0.4"
if (params.version) exit 0, log.info("$version")
params.download_db = false
params.skip_primer_trim = false
params.skip_nb = false
params.run_picrust2 = false
params.filterQ = 20
params.min_len = 1000
params.max_len = 1600
params.colorby = "condition"
params.skip_phylotree = false
params.minQ = 0
// Check input
params.input = false
if (params.input){
  n_sample=file(params.input).countLines() - 1
  if (n_sample == 1) {
    params.min_asv_totalfreq = 0
    params.min_asv_sample = 0
    println("Only 1 sample. min_asv_sample and min_asv_totalfreq set to 0.")
  } else {
    params.min_asv_totalfreq = 5
    params.min_asv_sample = 1
  }
} else {
  println("No input file given to --input!")
  n_sample=0
  params.min_asv_totalfreq = 0
  params.min_asv_sample = 0
}

if (params.skip_primer_trim) {
  params.front_p = 'none'
  params.adapter_p = 'none'
  trim_cutadapt = "No"
} else {
  trim_cutadapt = "Yes"
  // These are default V1-V9 adapter
  params.front_p = 'AGRGTTYGATYMTGGCTCAG'
  params.adapter_p = 'AAGTCGTAACAAGGTARCY'
}
// Hidden parameters in case need to trim with DADA2 removePrimers function
// StrainID primers
// params.front_p = 'AGRRTTYGATYHTDGYTYAG'
// params.adapter_p = 'AGTACYRHRARGGAANGR'
params.pooling_method = 'pseudo'
params.vsearch_db = "$projectDir/databases/GTDB_ssu_all_r207.qza"
params.vsearch_tax = "$projectDir/databases/GTDB_ssu_all_r207.taxonomy.qza"
params.silva_db = "$projectDir/databases/silva_nr99_v138.1_wSpecies_train_set.fa.gz"
params.refseq_db = "$projectDir/databases/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz"
params.gtdb_db = "$projectDir/databases/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz"
params.maxreject = 100
params.maxaccept = 100
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
params.publish_dir_mode = "symlink"

log_text = """
  Parameters set for pb-16S-nf pipeline for PacBio HiFi 16S
  =========================================================
  Number of samples in samples TSV: $n_sample
  Filter input reads above Q: $params.filterQ
  Trim primers with cutadapt: $trim_cutadapt
  Forward primer: $params.front_p
  Reverse primer: $params.adapter_p
  Minimum amplicon length filtered in DADA2: $params.min_len
  Maximum amplicon length filtered in DADA2: $params.max_len
  maxEE parameter for DADA2 filterAndTrim: $params.max_ee
  minQ parameter for DADA2 filterAndTrim: $params.minQ
  Pooling method for DADA2 denoise process: $params.pooling_method
  Minimum number of samples required to keep any ASV: $params.min_asv_sample
  Minimum number of reads required to keep any ASV: $params.min_asv_totalfreq 
  Taxonomy sequence database for VSEARCH: $params.vsearch_db
  Taxonomy annotation database for VSEARCH: $params.vsearch_tax
  Skip Naive Bayes classification: $params.skip_nb
  SILVA database for Naive Bayes classifier: $params.silva_db
  GTDB database for Naive Bayes classifier: $params.gtdb_db
  RefSeq + RDP database for Naive Bayes classifier: $params.refseq_db
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

log.info(log_text)

// Save parameters into file
process write_log{
  conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
  container "kpinpb/pb-16s-nf-tools:latest" 
  publishDir "$params.outdir", mode: params.publish_dir_mode

  input:
  val(logs)

  output:
  path "parameters.txt"

  script:
  """
  echo '$logs' > parameters.txt
  """
}

// QC before cutadapt
process QC_fastq {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
  container "kpinpb/pb-16s-nf-tools:latest"
  label 'cpu8'
  publishDir "$params.outdir/filtered_input_FASTQ", pattern: '*filterQ*.fastq.gz', mode: params.publish_dir_mode

  input:
  tuple val(sampleID), path(sampleFASTQ)

  output:
  path "${sampleID}.seqkit.readstats.tsv", emit: all_seqkit_stats
  path "${sampleID}.seqkit.summarystats.tsv", emit: all_seqkit_summary
  tuple val(sampleID), path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq
  path("${sampleID}.filterQ${params.filterQ}.fastq.gz"), emit: filtered_fastq_files

  script:
  """
  seqkit fx2tab -j $task.cpus -q --gc -l -H -n -i $sampleFASTQ |\
    csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.readstats.tsv
  seqkit stats -T -j $task.cpus -a ${sampleFASTQ} |\
    csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.summarystats.tsv
  seqkit seq -j $task.cpus --min-qual $params.filterQ $sampleFASTQ --out-file ${sampleID}.filterQ${params.filterQ}.fastq.gz
  """
}

// Trim full length 16S primers with cutadapt
process cutadapt {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/trimmed_primers_FASTQ", pattern: '*.fastq.gz', mode: params.publish_dir_mode
  publishDir "$params.outdir/cutadapt_summary", pattern: '*.report', mode: params.publish_dir_mode
  cpus params.cutadapt_cpu

  input:
  tuple val(sampleID), path(sampleFASTQ)

  output:
  tuple val(sampleID), path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq
  path "*.report", emit: cutadapt_summary
  path "cutadapt_summary_${sampleID}.tsv", emit: summary_tocollect
  path("${sampleID}.trimmed.fastq.gz"), emit: cutadapt_fastq_files

  script:
  """
  cutadapt -g "${params.front_p}...${params.adapter_p}" \
    ${sampleFASTQ} \
    -o ${sampleID}.trimmed.fastq.gz \
    -j ${task.cpus} --trimmed-only --revcomp -e 0.1 \
    --json ${sampleID}.cutadapt.report

  input_read=`jq -r '.read_counts | .input' ${sampleID}.cutadapt.report`
  demux_read=`jq -r '.read_counts | .output' ${sampleID}.cutadapt.report`
  echo -e "sample\tinput_reads\tdemuxed_reads" > cutadapt_summary_${sampleID}.tsv
  echo -e "${sampleID}\t\$input_read\t\$demux_read" >> cutadapt_summary_${sampleID}.tsv
  """
}

// QC after cutadapt
process QC_fastq_post_trim {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
  container "kpinpb/pb-16s-nf-tools:latest"
  label 'cpu8'
  publishDir "$params.outdir/filtered_input_FASTQ", pattern: '*post_trim.tsv', mode: params.publish_dir_mode

  input:
  tuple val(sampleID), path(sampleFASTQ)

  output:
  path "${sampleID}.seqkit.readstats.post_trim.tsv", emit: all_seqkit_stats

  script:
  """
  seqkit fx2tab -j $task.cpus -q --gc -l -H -n -i $sampleFASTQ |\
    csvtk mutate2 -C '%' -t -n sample -e '"${sampleID}"' > ${sampleID}.seqkit.readstats.post_trim.tsv
  """
}

// Collect QC into single files
process collect_QC {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
  container "kpinpb/pb-16s-nf-tools:latest"
  publishDir "$params.outdir/results/reads_QC", mode: params.publish_dir_mode
  label 'cpu8'

  input:
  path "*"
  path "*"
  path "*"
  path "*"

  output:
  path "all_samples_seqkit.readstats.tsv", emit: all_samples_readstats
  path "all_samples_seqkit.readstats.post_trim.tsv", emit: all_samples_readstats_post_trim
  path "all_samples_seqkit.summarystats.tsv", emit: all_samples_summarystats
  path "seqkit.summarised_stats.group_by_samples.tsv", emit: summarised_sample_readstats
  path "seqkit.summarised_stats.group_by_samples.pretty.tsv"
  path "all_samples_cutadapt_stats.tsv", emit: cutadapt_summary

  script:
  """
  csvtk concat -t -C '%' *.seqkit.readstats.tsv > all_samples_seqkit.readstats.tsv
  csvtk concat -t -C '%' *.seqkit.readstats.post_trim.tsv > all_samples_seqkit.readstats.post_trim.tsv
  csvtk concat -t -C '%' *.seqkit.summarystats.tsv > all_samples_seqkit.summarystats.tsv
  csvtk concat -t cutadapt_summary*.tsv > all_samples_cutadapt_stats.tsv
  # Summary read_qual for each sample
  csvtk summary -t -C '%' -g sample -f length:q1,length:q3,length:median,avg.qual:q1,avg.qual:q3,avg.qual:median all_samples_seqkit.readstats.tsv > seqkit.summarised_stats.group_by_samples.tsv
  csvtk pretty -t -C '%' seqkit.summarised_stats.group_by_samples.tsv > seqkit.summarised_stats.group_by_samples.pretty.tsv
  """
}

// Collect QC into single files if skipping cutadapt
process collect_QC_skip_cutadapt {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
  container "kpinpb/pb-16s-nf-tools:latest"
  publishDir "$params.outdir/results/reads_QC", mode: params.publish_dir_mode
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
  publishDir "$params.outdir/results/", mode: params.publish_dir_mode
  container "kpinpb/pb-16s-nf-tools:latest"

  input: 
  val(samplesFASTQ)
  path(metadata)

  output:
  path "samplefile*.txt", emit: sample_trimmed_file

  """
  echo -e "sample-id\tabsolute-filepath" > samplefile.txt
  echo -e "$samplesFASTQ" | sed -e 's/\\[//g' -e 's/\\]//g' | sed -r 's/\s+//g' | tr ',' '\\n' | split -l 2 - sample_ind_
  for i in \$(ls sample_ind_*); do sample=\$(cat \${i} | tr '\\n' '\\t' | cut -f1); fastq=\$(basename \$(cat \${i} | tr '\\n' '\\t' | cut -f2)); echo -e "\${sample}\t\${fastq}"; done >> samplefile.txt
  rm -f sample_ind_*
  # If pool column exists, split sample files
  poolyes=\$(csvtk headers -t $metadata | grep pool | wc -l)
  if [[ \$poolyes -eq 1 ]]
  then
    csvtk join -t samplefile.txt <(csvtk cut -t -f sample_name,pool $metadata) -f "sample-id;sample_name" | csvtk split -t -f pool
    counter=1
    for i in \$(ls stdin*.tsv); do csvtk cut -t -f -pool \${i} > samplefile\${counter}.txt; rm -f \${i}; let counter++; done
    rm -f samplefile.txt
  fi
  """

}

// Put all trimmed samples into a single file for QIIME2 import if skipping cutadapt
process prepare_qiime2_manifest_skip_cutadapt {
  label 'cpu_def'
  publishDir "$params.outdir/results/", mode: params.publish_dir_mode
  container "kpinpb/pb-16s-nf-tools:latest"

  input: 
  val(samplesFASTQ)
  path(metadata)

  output:
  path "samplefile.txt", emit: sample_trimmed_file

  """
  echo -e "sample-id\tabsolute-filepath" > samplefile.txt
  echo -e "$samplesFASTQ" | sed -e 's/\\[//g' -e 's/\\]//g' | sed -r 's/\s+//g' | tr ',' '\\n' | split -l 2 - sample_ind_
  for i in \$(ls sample_ind_*); do sample=\$(cat \${i} | tr '\\n' '\\t' | cut -f1); fastq=\$(basename \$(cat \${i} | tr '\\n' '\\t' | cut -f2)); echo -e "\${sample}\t\${fastq}"; done >> samplefile.txt
  rm -f sample_ind_*
  # If pool column exists, split sample files
  poolyes=\$(csvtk headers -t $metadata | grep pool | wc -l)
  if [[ \$poolyes -eq 1 ]]
  then
    csvtk join -t samplefile.txt <(csvtk cut -t -f sample_name,pool $metadata) -f "sample-id;sample_name" | csvtk split -t -f pool
    counter=1
    for i in \$(ls stdin*.tsv); do csvtk cut -t -f -pool \${i} > samplefile\${counter}.txt; rm -f \${i}; let counter++; done
    rm -f samplefile.txt
  fi
  """

}

// Import data into QIIME 2, note the `path *` in input is staging all the input 
// FASTQ into the process so it'll work on AWS that has no notion of local directory path
process import_qiime2 {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/import_qiime", mode: params.publish_dir_mode
  label 'cpu_def'

  input:
  path sample_manifest
  path '*'

  output:
  path 'samples.qza'

  script:
  """
  # Make sure path is absolute 
  awk -v wdir="\$(pwd)/" -F\$'\t' '{if (NR>1){print \$1,wdir\$2} else {print \$0}}' OFS=\$'\t' $sample_manifest > sample_list.txt
  qiime tools import --type 'SampleData[SequencesWithQuality]' \
    --input-path sample_list.txt \
    --output-path samples.qza \
    --input-format SingleEndFastqManifestPhred33V2
  """
}

// Merge sample manifest file for final vis purpose later on
process merge_sample_manifest {
  label 'cpu_def'
  container "kpinpb/pb-16s-nf-tools:latest"

  input:
  path sample_manifest

  output:
  path 'merged_sample_file.txt', emit: merged_samplefile

  script:
  """
  csvtk concat -t $sample_manifest > merged_sample_file.txt
  """
}

// Summarize demultiplex statistics
process demux_summarize {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/summary_demux", mode: params.publish_dir_mode
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

// Denoise into ASVs
process dada2_denoise {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
  cpus params.dada2_cpu

  input:
  path samples_qza
  path dada_ccs_script
  val(minQ)

  output:
  path "dada2-ccs_rep*.qza", emit: asv_seq
  path "dada2-ccs_table*.qza", emit: asv_freq
  path "dada2-ccs_stats*.qza", emit:asv_stats
  path "seqtab_nochim*.rds", emit: dada2_rds
  path "plot_error_model*.pdf"

  script:
  """
  # Use custom script that can skip primer trimming
  mkdir -p dada2_custom_script
  cp $dada_ccs_script dada2_custom_script/run_dada_ccs_original.R
  sed 's/minQ\\ =\\ 0/minQ=$minQ/g' dada2_custom_script/run_dada_ccs_original.R > \
    dada2_custom_script/run_dada_ccs.R
  chmod +x dada2_custom_script/run_dada_ccs.R
  export PATH="./dada2_custom_script:\$PATH"
  which run_dada_ccs.R
  qiime dada2 denoise-ccs --i-demultiplexed-seqs $samples_qza \
    --o-table dada2-ccs_table.qza \
    --o-representative-sequences dada2-ccs_rep.qza \
    --o-denoising-stats dada2-ccs_stats.qza \
    --p-min-len $params.min_len --p-max-len $params.max_len \
    --p-max-ee $params.max_ee \
    --p-front 'none' \
    --p-adapter 'none' \
    --p-n-threads $task.cpus \
    --p-pooling-method \'$params.pooling_method\'

  """
}

// TODO Merge ASV table and sequences of all chunks, if chunked
// Take note of the output above! There will only be one set of out from
// different group in dada2 folder after dada2_denoise
process mergeASV {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
  cpus params.dada2_cpu

  input:
  path "dada2-ccs_rep*.qza"
  path "dada2-ccs_table*.qza"
  path "dada2-ccs_stats*.qza"

  output:
  path "dada2-ccs_rep_merged.qza", emit: asv_seq
  path "dada2-ccs_table_merged.qza", emit: asv_freq
  path "dada2-ccs_stats_merged.qza", emit: asv_stats

  script:
  """
  qiime feature-table merge \
    --i-tables dada2-ccs_table*.qza \
    --o-merged-table dada2-ccs_table_merged.qza
  
  qiime feature-table merge-seqs \
    --i-data dada2-ccs_rep*.qza \
    --o-merged-data dada2-ccs_rep_merged.qza

  for i in \$(ls dada2-ccs_stats*.qza);
  do
    qiime tools export --input-path \${i} \
    --output-path ./\${i%%.qza}_export/
  done

  cat ./*_export/stats.tsv | awk '{if (NR>2 && \$0~/^#|sample-id/) {next} else {print \$0}}' > merged_stats.tsv
  qiime tools import --input-path merged_stats.tsv --output-path dada2-ccs_stats_merged --type 'SampleData[DADA2Stats]'
  """
}


// Filter DADA2 ASVs using minimum number of samples and frequency of ASVs
process filter_dada2 {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/dada2", mode: params.publish_dir_mode
  cpus params.dada2_cpu

  input:
  path asv_table
  path asv_seq

  output:
  path "dada2-ccs_rep_filtered.qza", emit: asv_seq
  path "dada2-ccs_table_filtered.qza", emit: asv_freq
  path "dada2_ASV.fasta", emit:asv_seq_fasta

  script:
  """
  if [ $params.min_asv_sample -gt 0 ] && [ $params.min_asv_totalfreq -gt 0 ]
  then
    # Filter ASVs and sequences
    qiime feature-table filter-features \
      --i-table $asv_table \
      --p-min-frequency $params.min_asv_totalfreq \
      --p-min-samples $params.min_asv_sample \
      --p-filter-empty-samples \
      --o-filtered-table dada2-ccs_table_filtered.qza
  elif [ $params.min_asv_sample -gt 0 ] && [ $params.min_asv_totalfreq -eq 0 ]
  then
    qiime feature-table filter-features \
      --i-table $asv_table \
      --p-min-samples $params.min_asv_sample \
      --p-filter-empty-samples \
      --o-filtered-table dada2-ccs_table_filtered.qza
  elif [ $params.min_asv_sample -eq 0 ] && [ $params.min_asv_totalfreq -gt 0 ]
  then
    qiime feature-table filter-features \
      --i-table $asv_table \
      --p-min-frequency $params.min_asv_totalfreq \
      --p-filter-empty-samples \
      --o-filtered-table dada2-ccs_table_filtered.qza
  else
    mv $asv_table dada2-ccs_table_filtered.qza
  fi
  qiime feature-table filter-seqs \
    --i-data $asv_seq \
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
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", pattern: 'best_tax*', mode: params.publish_dir_mode
  publishDir "$params.outdir/nb_tax", mode: params.publish_dir_mode
  cpus params.vsearch_cpu

  input:
  path asv_seq_fasta
  path asv_seq
  path asv_freq
  path silva_db
  path gtdb_db
  path refseq_db
  path assign_script

  output:
  path "best_tax.qza", emit:best_nb_tax_qza
  path "best_taxonomy.tsv", emit: best_nb_tax
  path "best_taxonomy_withDB.tsv"
  path "best_tax_merged_freq_tax.tsv", emit: best_nb_tax_tsv
  path "silva_nb.tsv"
  path "gtdb_nb.tsv"
  path "refseq_rdp_nb.tsv"

  script:
  """
  Rscript --vanilla $assign_script $asv_seq_fasta $task.cpus $silva_db $gtdb_db $refseq_db 80

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
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
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
  env(int_alpha_d), emit:alpha_depth

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

  # Get number of reads for ASV covering 80% of samples
  number=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | wc -l`
  ninety=0.8
  # Handle samples less than 5 as we are using at least 20% of samples for
  # rarefaction sampling depth (i.e. if sample is less than 5 the math result is 0)
  if [ \${number} -le 2 ];
  then
    rarefaction_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
      tail -n1 | cut -f2 -d,`
    int_rarefaction_d=\${rarefaction_d%%.*}
    echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt
    alpha_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
      tail -n1 | cut -f2 -d,`
    int_alpha_d=\${alpha_d%%.*}
  elif [ \${number} -gt 2 ] && [ \${number} -lt 5 ];
  then
    result=`echo "(\$number * \$ninety)" | bc -l` 
    rarefaction_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
      head -n \${result%%.*} | tail -n1 | cut -f2 -d,`
    int_rarefaction_d=\${rarefaction_d%%.*}
    echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt
    alpha_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
      tail -n1 | cut -f2 -d,`
    int_alpha_d=\${alpha_d%%.*}
  else
    result=`echo "(\$number * \$ninety)" | bc -l` 
    rarefaction_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
      head -n \${result%%.*} | tail -n1 | cut -f2 -d,`
    int_rarefaction_d=\${rarefaction_d%%.*}
    echo \${int_rarefaction_d} > rarefaction_depth_suggested.txt
    # Use lowest depth covering just 20% of sample for alpha rarefaction curve
    alpha_res=`echo "(\$number * 0.2)" | bc -l`
    alpha_d=`tail -n+2 dada2_table_summary/sample-frequency-detail.csv | sort -t, -k2 -nr | \
      head -n \${alpha_res%%.*} | tail -n1 | cut -f2 -d,`
    int_alpha_d=\${alpha_d%%.*}
  fi

  qiime tools export --input-path dada2_stats.qzv \
    --output-path ./dada2_stats
  mv dada2_stats/metadata.tsv dada2_qc.tsv
  """
}

process qiime2_phylogeny_diversity {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results/phylogeny_diversity", mode: params.publish_dir_mode
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
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
  label 'cpu_def'

  input:
  path asv_freq
  path metadata
  val(alpha_d)

  output:
  path "*"

  script:
  """
  qiime diversity alpha-rarefaction --i-table $asv_freq \
    --m-metadata-file $metadata \
    --o-visualization alpha-rarefaction-curves.qzv \
    --p-min-depth 10 --p-max-depth $alpha_d
  """
}

// Classify taxonomy and export table using VSEARCH approach
process class_tax {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
  cpus params.vsearch_cpu

  input:
  path asv_seq
  path asv_freq
  path vsearch_db
  path vsearch_tax

  output:
  path "taxonomy.vsearch.qza", emit: tax_vsearch
  path "tax_export/taxonomy.tsv", emit: tax_tsv
  path "merged_freq_tax.qzv", emit: tax_freq_tab
  path "vsearch_merged_freq_tax.tsv", emit: tax_freq_tab_tsv

  script:
  """
  qiime feature-classifier classify-consensus-vsearch --i-query $asv_seq \
    --o-classification taxonomy.vsearch.qza \
    --i-reference-reads $vsearch_db \
    --i-reference-taxonomy $vsearch_tax \
    --p-threads $task.cpus \
    --p-maxrejects $params.maxreject \
    --p-maxaccepts $params.maxaccept \
    --p-perc-identity $params.vsearch_identity \
    --p-top-hits-only

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
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
  label 'cpu_def'

  input:
  path asv_freq
  path tax_tsv
  path tax_tsv_vsearch

  output:
  path "feature-table-tax.biom", emit:biom
  path "feature-table-tax_vsearch.biom", emit:biom_vsearch

  script:
  """
  qiime tools export --input-path $asv_freq --output-path asv_freq/

  sed 's/Feature ID/#OTUID/' $tax_tsv | sed 's/Taxon/taxonomy/' | \
    sed 's/Consensus/confidence/' > biom-taxonomy.tsv

  biom add-metadata -i asv_freq/feature-table.biom \
    -o feature-table-tax.biom \
    --observation-metadata-fp biom-taxonomy.tsv \
    --sc-separated taxonomy

  sed 's/Feature ID/#OTUID/' $tax_tsv_vsearch | sed 's/Taxon/taxonomy/' | \
    sed 's/Consensus/confidence/' > biom-taxonomy_vsearch.tsv

  biom add-metadata -i asv_freq/feature-table.biom \
    -o feature-table-tax_vsearch.biom \
    --observation-metadata-fp biom-taxonomy_vsearch.tsv \
    --sc-separated taxonomy
  """
}

// Export results into biom for use with phyloseq
process export_biom_skip_nb {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
  label 'cpu_def'

  input:
  path asv_freq
  path tax_tsv_vsearch

  output:
  path "feature-table-tax_vsearch.biom", emit:biom_vsearch

  script:
  """
  qiime tools export --input-path $asv_freq --output-path asv_freq/

  sed 's/Feature ID/#OTUID/' $tax_tsv_vsearch | sed 's/Taxon/taxonomy/' | \
    sed 's/Consensus/confidence/' > biom-taxonomy_vsearch.tsv

  biom add-metadata -i asv_freq/feature-table.biom \
    -o feature-table-tax_vsearch.biom \
    --observation-metadata-fp biom-taxonomy_vsearch.tsv \
    --sc-separated taxonomy
  """
}

process picrust2 {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-pbtools.yml" : null)
  container "kpinpb/pb-16s-nf-tools:latest"
  label 'cpu32'
  publishDir "$params.outdir/results", mode: params.publish_dir_mode

  input:
  path dada2_asv
  path vsearch_biom

  output:
  path "picrust2"

  script:
  """
  picrust2_pipeline.py -s $dada2_asv -i $vsearch_biom \
    -o picrust2 \
    --in_traits EC,KO \
    --verbose -p $task.cpus
  """
}

process barplot {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
  label 'cpu_def'

  input:
  path asv_tab
  path tax
  path metadata
  val(oname)

  output:
  path "${oname}"

  script:
  """
  qiime taxa barplot --i-table $asv_tab --i-taxonomy $tax \
    --m-metadata-file $metadata \
    --o-visualization $oname
  """
}

process barplot_nb {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
  label 'cpu_def'

  input:
  path asv_tab
  path tax
  path metadata
  val(oname)

  output:
  path "${oname}"

  script:
  """
  qiime taxa barplot --i-table $asv_tab --i-taxonomy $tax \
    --m-metadata-file $metadata \
    --o-visualization $oname
  """
}

// HTML report
process html_rep {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-vis-conda.yml" : null)
  container "kpinpb/pb-16s-vis:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
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
  path post_trim_readstats
  path rmd_vis_biom_script
  path rmd_helper

  output:
  path "visualize_biom.html", emit: html_report

  script:
  if (params.enable_container)
  """
  export R_LIBS_USER="/opt/conda/envs/pb-16S-vis/lib/R/library"
  # Copy script because of path issue (Rmd expects input files to be in the same folder as Rmd file)
  cp $rmd_vis_biom_script vis_script.rmd
  Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
  mv vis_script.html visualize_biom.html
  """
  else
  """
  cp $rmd_vis_biom_script vis_script.rmd
  # Copy script because of path issue (Rmd expects input files to be in the same folder as Rmd file)
  Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
  mv vis_script.html visualize_biom.html
  """
}

// HTML report
process html_rep_skip_cutadapt {
  conda (params.enable_conda ? "$projectDir/env/pb-16s-vis-conda.yml" : null)
  container "kpinpb/pb-16s-vis:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
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
  val(post_trim_readstats)
  path rmd_vis_biom_script
  path rmd_helper

  output:
  path "visualize_biom.html", emit: html_report

  script:
  if (params.enable_container)
  """
  export R_LIBS_USER="/opt/conda/envs/pb-16S-vis/lib/R/library"
  # Copy script because of path issue (Rmd expects input files to be in the same folder as Rmd file)
  cp $rmd_vis_biom_script vis_script.rmd
  Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
  mv vis_script.html visualize_biom.html
  """
  else
  """
  # Copy script because of path issue (Rmd expects input files to be in the same folder as Rmd file)
  cp $rmd_vis_biom_script vis_script.rmd
  Rscript -e 'rmarkdown::render("vis_script.rmd", params=list(merged_tax_tab_file="$tax_freq_tab_tsv", metadata="$metadata", sample_file="$sample_manifest", dada2_qc="$dada2_qc", reads_qc="$reads_qc", summarised_reads_qc="$summarised_reads_qc", cutadapt_qc="$cutadapt_summary_qc", vsearch_tax_tab_file="$vsearch_tax_tsv", colorby="$colorby", bray_mat="$bray_mat", unifrac_mat="$unifrac_mat", wunifrac_mat="$wunifrac_mat", post_trim_readstats="$post_trim_readstats"), output_dir="./")'
  mv vis_script.html visualize_biom.html
  """
}

// Krona plot
process krona_plot {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$params.outdir/results", mode: params.publish_dir_mode
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
  if (params.enable_conda){
  """
    pip install git+https://github.com/kaanb93/q2-krona.git
    qiime krona collapse-and-plot --i-table $asv_freq --i-taxonomy $taxonomy --o-krona-plot krona.qzv
    qiime tools export --input-path krona.qzv --output-path krona_html
  """
  } else {
  """
    qiime krona collapse-and-plot --i-table $asv_freq --i-taxonomy $taxonomy --o-krona-plot krona.qzv
    qiime tools export --input-path krona.qzv --output-path krona_html
  """
  }
}

process download_db {
  conda (params.enable_conda ? "$projectDir/env/qiime2-2022.2-py38-linux-conda.yml" : null)
  container "kpinpb/pb-16s-nf-qiime:latest"
  publishDir "$projectDir/databases", mode: "copy"
  label 'cpu_def'

  output:
  path "*"

  script:
  """
  echo "Downloading SILVA sequences for VSEARCH..."
  wget -N --content-disposition 'https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1'
  echo "Downloading GTDB sequences for VSEARCH..."
  wget -N --content-disposition 'https://zenodo.org/record/6655692/files/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz?download=1'
  echo "Downloading RefSeq + RDP sequences for VSEARCH..."
  wget -N --content-disposition 'https://zenodo.org/record/4735821/files/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz?download=1'

  echo "Downloading SILVA sequences and taxonomies for Naive Bayes classifier"
  wget -N --content-disposition 'https://data.qiime2.org/2022.2/common/silva-138-99-seqs.qza'
  wget -N --content-disposition 'https://data.qiime2.org/2022.2/common/silva-138-99-tax.qza'

  echo "Downloading GTDB database"
  wget -N --content-disposition --no-check-certificate 'https://zenodo.org/record/6912512/files/GTDB_ssu_all_r207.taxonomy.qza?download=1'
  wget -N --content-disposition --no-check-certificate 'https://zenodo.org/record/6912512/files/GTDB_ssu_all_r207.qza?download=1'
  #   wget -N --content-disposition --no-check-certificate 'https://data.gtdb.ecogenomic.org/releases/release207/207.0/genomic_files_all/ssu_all_r207.tar.gz'
  #   tar xfz ssu_all_r207.tar.gz
  #   mv ssu_all_r207.fna GTDB_ssu_all_r207.fna
  #   grep "^>" GTDB_ssu_all_r207.fna | awk '{print gensub(/^>(.*)/, "\\\\1", "g", \$1),gensub(/^>.*\\ (d__.*)\\ \\[.*\\[.*\\[.*/, "\\\\1", "g", \$0)}' OFS=\$'\\t' > GTDB_ssu_all_r207.taxonomy.tsv
  #   qiime tools import --type 'FeatureData[Sequence]' --input-path GTDB_ssu_all_r207.fna --output-path GTDB_ssu_all_r207.qza
  #   qiime tools import --type 'FeatureData[Taxonomy]' --input-format HeaderlessTSVTaxonomyFormat \
  #     --input-path GTDB_ssu_all_r207.taxonomy.tsv --output-path GTDB_ssu_all_r207.taxonomy.qza
  """
}

workflow pb16S {
  if (params.download_db){
    download_db()
  } else {
    write_log(log_text)
    sample_file = channel.fromPath(params.input)
      .splitCsv(header: ['sample', 'fastq'], skip: 1, sep: "\t")
      .map{ row -> tuple(row.sample, file(row.fastq)) }
    metadata_file = channel.fromPath(params.metadata)
    QC_fastq(sample_file)
    if (params.skip_primer_trim){
      collect_QC_skip_cutadapt(QC_fastq.out.all_seqkit_stats.collect(),
          QC_fastq.out.all_seqkit_summary.collect())
      prepare_qiime2_manifest_skip_cutadapt(QC_fastq.out.filtered_fastq.collect(), metadata_file)
      cutadapt_summary = "none"
      qiime2_manifest = prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file.flatten()
      merge_sample_manifest(prepare_qiime2_manifest_skip_cutadapt.out.sample_trimmed_file.collect())
      collect_QC_readstats = collect_QC_skip_cutadapt.out.all_samples_readstats
      post_trim_readstats = "none"
      collect_QC_summarised_sample_stats = collect_QC_skip_cutadapt.out.summarised_sample_readstats
      import_qiime2(qiime2_manifest, QC_fastq.out.filtered_fastq_files.collect())
    } else {
      cutadapt(QC_fastq.out.filtered_fastq)
      QC_fastq_post_trim(cutadapt.out.cutadapt_fastq)
      collect_QC(QC_fastq.out.all_seqkit_stats.collect(),
          QC_fastq.out.all_seqkit_summary.collect(),
          cutadapt.out.summary_tocollect.collect(),
          QC_fastq_post_trim.out.all_seqkit_stats.collect())
      post_trim_readstats = collect_QC.out.all_samples_readstats_post_trim
      prepare_qiime2_manifest(cutadapt.out.cutadapt_fastq.collect(), metadata_file)
      cutadapt_summary = collect_QC.out.cutadapt_summary
      qiime2_manifest = prepare_qiime2_manifest.out.sample_trimmed_file.flatten()
      merge_sample_manifest(prepare_qiime2_manifest.out.sample_trimmed_file.collect())
      collect_QC_readstats = collect_QC.out.all_samples_readstats
      collect_QC_summarised_sample_stats = collect_QC.out.summarised_sample_readstats
      import_qiime2(qiime2_manifest, cutadapt.out.cutadapt_fastq_files.collect())
    }
    demux_summarize(import_qiime2.out)
    dada2_denoise(import_qiime2.out, params.dadaCCS_script, params.minQ)
    mergeASV(dada2_denoise.out.asv_seq.collect(), dada2_denoise.out.asv_freq.collect(), dada2_denoise.out.asv_stats.collect())
    filter_dada2(mergeASV.out.asv_freq, mergeASV.out.asv_seq)
    dada2_qc(mergeASV.out.asv_stats, filter_dada2.out.asv_freq, metadata_file)
    if( params.rarefaction_depth > 0 ){
      rd = params.rarefaction_depth
      qiime2_phylogeny_diversity(metadata_file, filter_dada2.out.asv_seq,
          filter_dada2.out.asv_freq, rd)
    } else {
      qiime2_phylogeny_diversity(metadata_file, filter_dada2.out.asv_seq,
          filter_dada2.out.asv_freq, dada2_qc.out.rarefaction_depth)
    }
    dada2_rarefaction(filter_dada2.out.asv_freq, metadata_file, dada2_qc.out.alpha_depth)
    class_tax(filter_dada2.out.asv_seq, filter_dada2.out.asv_freq, params.vsearch_db, params.vsearch_tax)
    if(params.skip_nb){
      nb_tax = []
      export_biom_skip_nb(filter_dada2.out.asv_freq, class_tax.out.tax_tsv) 
      barplot(filter_dada2.out.asv_freq, class_tax.out.tax_vsearch, metadata_file, "taxonomy_barplot_vsearch.qzv")
    } else {
      dada2_assignTax(filter_dada2.out.asv_seq_fasta, filter_dada2.out.asv_seq, filter_dada2.out.asv_freq,
          params.silva_db, params.gtdb_db, params.refseq_db, params.dadaAssign_script)
      nb_tax = dada2_assignTax.out.best_nb_tax_tsv 
      export_biom(filter_dada2.out.asv_freq, dada2_assignTax.out.best_nb_tax, class_tax.out.tax_tsv)
      barplot_nb(filter_dada2.out.asv_freq, dada2_assignTax.out.best_nb_tax_qza, metadata_file, "taxonomy_barplot_nb.qzv")
      barplot(filter_dada2.out.asv_freq, class_tax.out.tax_vsearch, metadata_file, "taxonomy_barplot_vsearch.qzv")
    }
    if (params.run_picrust2){
      picrust2(filter_dada2.out.asv_seq_fasta, export_biom.out.biom_vsearch)
    }
    if (params.skip_primer_trim){
      html_rep_skip_cutadapt(nb_tax, metadata_file, merge_sample_manifest.out.merged_samplefile,
          dada2_qc.out.dada2_qc_tsv, 
          collect_QC_readstats, collect_QC_summarised_sample_stats,
          cutadapt_summary, class_tax.out.tax_freq_tab_tsv, qiime2_phylogeny_diversity.out.bray_mat,
          qiime2_phylogeny_diversity.out.unifrac_mat, qiime2_phylogeny_diversity.out.wunifrac_mat,
          params.colorby, post_trim_readstats, params.rmd_vis_biom_script, params.rmd_helper)
    } else {
      html_rep(nb_tax, metadata_file,  merge_sample_manifest.out.merged_samplefile,
          dada2_qc.out.dada2_qc_tsv, 
          collect_QC_readstats, collect_QC_summarised_sample_stats,
          cutadapt_summary, class_tax.out.tax_freq_tab_tsv, qiime2_phylogeny_diversity.out.bray_mat,
          qiime2_phylogeny_diversity.out.unifrac_mat, qiime2_phylogeny_diversity.out.wunifrac_mat,
          params.colorby, post_trim_readstats, params.rmd_vis_biom_script, params.rmd_helper )
    }
    krona_plot(filter_dada2.out.asv_freq, class_tax.out.tax_vsearch)
  }
}

workflow {
  pb16S()
}
