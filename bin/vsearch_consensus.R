#!/usr/bin/env Rscript

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  err_quit(
    paste(
      "Expected 3 arguments:",
      "1) input_vsearch_hits_tsv",
      "2) output_consensus_tsv",
      "3) min_consensus",
      sep = "\n"
    )
  )
}

input_tsv <- args[1]
output_tsv <- args[2]
min_consensus <- suppressWarnings(as.numeric(args[3]))

if (!file.exists(input_tsv)) {
  err_quit(sprintf("Input file does not exist: %s", input_tsv))
}

if (is.na(min_consensus) || min_consensus < 0 || min_consensus > 1) {
  err_quit("min_consensus must be a number between 0 and 1")
}


# -------------------------------------------------------------------------
# Read VSEARCH top hits
# -------------------------------------------------------------------------

hits <- read.table(
  input_tsv,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c(
  "Feature ID",
  "Reference ID",
  "Identity",
  "AlignmentLength",
  "Taxon"
)

missing_cols <- setdiff(required_cols, colnames(hits))

if (length(missing_cols) > 0) {
  err_quit(
    sprintf(
      "Input file is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (nrow(hits) == 0) {
  write.table(
    data.frame(
      `Feature ID` = character(),
      Taxon = character(),
      Confidence = numeric(),
      Identity = numeric(),
      AlignmentLength = numeric(),
      TopHits = integer(),
      check.names = FALSE
    ),
    file = output_tsv,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  quit(save = "no", status = 0)
}


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

split_taxonomy <- function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  trimws(parts)
}


consensus_taxonomy <- function(taxa, min_consensus) {

  taxa_split <- lapply(taxa, split_taxonomy)

  max_rank <- max(lengths(taxa_split))

  consensus <- character()
  last_support <- NA_real_

  for (rank_idx in seq_len(max_rank)) {

    rank_values <- vapply(
      taxa_split,
      function(x) {
        if (length(x) >= rank_idx) {
          x[rank_idx]
        } else {
          NA_character_
        }
      },
      character(1)
    )

    rank_values <- rank_values[
      !is.na(rank_values) &
      nzchar(rank_values) &
      rank_values != "Unclassified"
    ]

    if (length(rank_values) == 0) {
      break
    }

    counts <- sort(
      table(rank_values),
      decreasing = TRUE
    )

    best_count <- as.integer(counts[1])
    best_taxon <- names(counts)[1]

    # Support is calculated against all retained top hits,
    # not only those with an annotation at this rank.
    support <- best_count / length(taxa)

    if (support < min_consensus) {
      break
    }

    consensus <- c(consensus, best_taxon)
    last_support <- support
  }

  if (length(consensus) == 0) {
    return(
      list(
        taxonomy = "Unclassified",
        confidence = 0
      )
    )
  }

  list(
    taxonomy = paste(consensus, collapse = ";"),
    confidence = last_support
  )
}


# -------------------------------------------------------------------------
# Calculate consensus per ASV
# -------------------------------------------------------------------------

feature_ids <- unique(hits[["Feature ID"]])

results <- lapply(
  feature_ids,
  function(feature_id) {

    x <- hits[hits[["Feature ID"]] == feature_id, , drop = FALSE]

    consensus <- consensus_taxonomy(
      x[["Taxon"]],
      min_consensus
    )

    data.frame(
      `Feature ID` = feature_id,
      Taxon = consensus$taxonomy,
      Confidence = consensus$confidence,
      Identity = max(as.numeric(x[["Identity"]]), na.rm = TRUE),
      AlignmentLength = max(as.numeric(x[["AlignmentLength"]]), na.rm = TRUE),
      TopHits = nrow(x),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
)

result_df <- do.call(
  rbind,
  results
)

rownames(result_df) <- NULL


# -------------------------------------------------------------------------
# Write result
# -------------------------------------------------------------------------

write.table(
  result_df,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

quit(save = "no", status = 0)
