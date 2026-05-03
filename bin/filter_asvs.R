#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

parse_required_numeric <- function(x, arg_name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value) || value < 0) {
    err_quit(sprintf("Argument '%s' must be a non-negative number. Got: %s", arg_name, x))
  }
  value
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  err_quit(
    paste(
      "Expected 7 positional arguments:",
      "1) input_seqtab_nochim_rds",
      "2) output_seqtab_filtered_rds",
      "3) output_table_tsv",
      "4) output_asv_fasta",
      "5) output_filter_stats_tsv",
      "6) min_asv_totalfreq",
      "7) min_asv_sample",
      sep = "\n"
    )
  )
}

input_seqtab_rds <- args[1]
output_seqtab_rds <- args[2]
output_table_tsv <- args[3]
output_asv_fasta <- args[4]
output_filter_stats_tsv <- args[5]
min_asv_totalfreq <- parse_required_numeric(args[6], "min_asv_totalfreq")
min_asv_sample <- parse_required_numeric(args[7], "min_asv_sample")

if (!file.exists(input_seqtab_rds)) {
  err_quit(sprintf("Input seqtab RDS does not exist: %s", input_seqtab_rds))
}

seqtab <- readRDS(input_seqtab_rds)

if (!is.matrix(seqtab)) {
  err_quit("Loaded seqtab is not a matrix.")
}

if (nrow(seqtab) == 0 || ncol(seqtab) == 0) {
  err_quit("Input seqtab has zero rows or zero columns.")
}

input_n_samples <- nrow(seqtab)
input_n_asvs <- ncol(seqtab)
input_total_reads <- sum(seqtab)

keep_totalfreq <- if (min_asv_totalfreq > 0) {
  colSums(seqtab) >= min_asv_totalfreq
} else {
  rep(TRUE, ncol(seqtab))
}

keep_samplecount <- if (min_asv_sample > 0) {
  colSums(seqtab > 0) >= min_asv_sample
} else {
  rep(TRUE, ncol(seqtab))
}

keep_asvs <- keep_totalfreq & keep_samplecount
seqtab_filt <- seqtab[, keep_asvs, drop = FALSE]

# Drop empty samples, matching old QIIME behavior
keep_samples <- rowSums(seqtab_filt) > 0
seqtab_filt <- seqtab_filt[keep_samples, , drop = FALSE]

if (ncol(seqtab_filt) == 0) {
  err_quit("All ASVs were removed by filtering.")
}

saveRDS(seqtab_filt, output_seqtab_rds)

write.table(
  seqtab_filt,
  file = output_table_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)
asv_seqs <- colnames(seqtab_filt)
repseqs <- DNAStringSet(asv_seqs)
names(repseqs) <- asv_seqs
writeXStringSet(repseqs, filepath = output_asv_fasta)


stats_df <- data.frame(
  input_samples = input_n_samples,
  output_samples = nrow(seqtab_filt),
  removed_samples = input_n_samples - nrow(seqtab_filt),
  input_asvs = input_n_asvs,
  output_asvs = ncol(seqtab_filt),
  removed_asvs = input_n_asvs - ncol(seqtab_filt),
  input_total_reads = input_total_reads,
  output_total_reads = sum(seqtab_filt),
  removed_reads = input_total_reads - sum(seqtab_filt),
  min_asv_totalfreq = min_asv_totalfreq,
  min_asv_sample = min_asv_sample,
  stringsAsFactors = FALSE
)

write.table(
  stats_df,
  file = output_filter_stats_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

quit(save = "no", status = 0)
