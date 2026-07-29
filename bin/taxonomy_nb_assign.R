#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  err_quit(
    paste(
      "Expected 4 positional arguments:",
      "1) asv_fasta",
      "2) db_fasta",
      "3) db_name",
      "4) n_threads",
      sep = "\n"
    )
  )
}

asv_fasta <- args[1]
db_fasta <- args[2]
db_name <- args[3]
n_threads <- suppressWarnings(as.integer(args[4]))

if (!file.exists(asv_fasta)) err_quit(sprintf("ASV FASTA not found: %s", asv_fasta))
if (!file.exists(db_fasta)) err_quit(sprintf("DB FASTA not found: %s", db_fasta))
if (is.na(n_threads) || n_threads < 1) n_threads <- 1L

seqs <- readDNAStringSet(asv_fasta)
seq_vec <- as.character(seqs)
feature_ids <- names(seqs)

if (length(feature_ids) == 0 || is.null(feature_ids)) {
  err_quit("No ASV identifiers found in FASTA headers")
}

message("Running assignTaxonomy for DB: ", db_name)

tax <- assignTaxonomy(
  seqs = seq_vec,
  refFasta = db_fasta,
  multithread = n_threads,
  tryRC = TRUE,
  outputBootstraps = TRUE,
  verbose = TRUE
)

taxa_mat <- tax$tax
boot_mat <- tax$boot

if (is.null(dim(taxa_mat))) {
  err_quit(sprintf("Unexpected assignTaxonomy output for DB: %s", db_name))
}

collapse_taxonomy <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return("Unclassified")
  }
  paste(x, collapse = ";")
}

confidence_fun <- function(taxa, boot) {
  assigned <- which(!is.na(taxa) & taxa != "")

  if (length(assigned) == 0) {
    return(NA_real_)
  }

  deepest_rank <- max(assigned)

  suppressWarnings(as.numeric(boot[deepest_rank]))
}

confidences <- vapply(
  seq_len(nrow(taxa_mat)),
  function(i) {
    confidence_fun(
      taxa_mat[i, ],
      boot_mat[i, ]
    )
  },
  numeric(1)
)

taxon_strings <- apply(taxa_mat, 1, collapse_taxonomy)

out <- data.frame(
  FeatureID = feature_ids,
  Taxon = taxon_strings,
  Confidence = confidences,
  stringsAsFactors = FALSE
)

out_file <- sprintf("%s_nb.tsv", db_name)

write.table(
  out,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("Wrote: ", out_file)
quit(save = "no", status = 0)
