#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: taxonomy_vsearch_best.R ",
    "<vsearch_counts_1.tsv> [vsearch_counts_2.tsv ...] <db_priority>"
  )
}

# Last argument is the database priority.
db_priority_arg <- args[length(args)]
tax_files <- args[-length(args)]

db_priority <- unlist(strsplit(db_priority_arg, ","))
db_priority <- trimws(db_priority)
db_priority <- db_priority[nzchar(db_priority)]

if (length(tax_files) == 0) {
  stop("No VSEARCH taxonomy files provided")
}

# ------------------------------------------------------------------------------
# Read and annotate all database-specific VSEARCH tables
# ------------------------------------------------------------------------------

tables <- lapply(tax_files, function(file) {

  dt <- fread(file, check.names = FALSE)

  required_cols <- c(
    "Sequence",
    "Taxon",
    "Confidence",
    "Identity",
    "AlignmentLength",
    "TopHits"
  )

  missing_cols <- setdiff(required_cols, colnames(dt))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "File '%s' is missing required columns: %s",
      file,
      paste(missing_cols, collapse = ", ")
    ))
  }

  db <- basename(file)
  db <- sub("_vsearch_counts\\.tsv$", "", db)

  dt[, Database := db]

  dt[, Sequence := as.character(Sequence)]
  dt[, Taxon := as.character(Taxon)]
  dt[, Confidence := as.numeric(Confidence)]
  dt[, Identity := as.numeric(Identity)]
  dt[, AlignmentLength := as.numeric(AlignmentLength)]
  dt[, TopHits := as.integer(TopHits)]

  dt
})

all_tax <- rbindlist(tables, use.names = TRUE, fill = TRUE)

# ------------------------------------------------------------------------------
# Database priority
# Lower number = higher priority
# ------------------------------------------------------------------------------

all_tax[, db_priority := match(Database, db_priority)]

# Databases not explicitly listed in the priority parameter come last.
all_tax[
  is.na(db_priority),
  db_priority := length(db_priority) + 1L
]

# ------------------------------------------------------------------------------
# Determine whether an ASV actually has a VSEARCH assignment
# ------------------------------------------------------------------------------

all_tax[
  ,
  has_hit := !is.na(Taxon) &
             Taxon != "" &
             Taxon != "Unclassified" &
             !is.na(Identity)
]

# ------------------------------------------------------------------------------
# Pick best assignment per ASV
#
# Priority:
#   1. VSEARCH hit over no hit
#   2. highest sequence identity
#   3. longest alignment
#   4. configured database priority
# ------------------------------------------------------------------------------

setorder(
  all_tax,
  Sequence,
  -has_hit,
  -Confidence,
  -Identity,
  -AlignmentLength,
  db_priority
)

best <- all_tax[, .SD[1], by = Sequence]

# ------------------------------------------------------------------------------
# Output with database information
# ------------------------------------------------------------------------------

best_with_db <- copy(best)

best_with_db[
  ,
  c("has_hit", "db_priority") := NULL
]

fwrite(
  best_with_db,
  "best_vsearch_taxonomy_withDB.tsv",
  sep = "\t",
  quote = FALSE,
  na = ""
)

# ------------------------------------------------------------------------------
# Output without database information
# ------------------------------------------------------------------------------

best_no_db <- copy(best_with_db)

best_no_db[, Database := NULL]

fwrite(
  best_no_db,
  "best_vsearch_taxonomy.tsv",
  sep = "\t",
  quote = FALSE,
  na = ""
)

message(
  sprintf(
    "Selected best VSEARCH assignments for %d ASVs across %d databases",
    nrow(best),
    length(tax_files)
  )
)
