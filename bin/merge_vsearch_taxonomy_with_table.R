#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: merge_vsearch_taxonomy_with_table.R <taxonomy.tsv> <asv_table.tsv> <out.tsv>")
}

tax_file <- args[1]
asv_file <- args[2]
out_file <- args[3]

if (!file.exists(tax_file)) {
  stop(sprintf("Taxonomy file not found: %s", tax_file))
}

if (!file.exists(asv_file)) {
  stop(sprintf("ASV table file not found: %s", asv_file))
}

# ------------------------------------------------------------------------------
# Read VSEARCH taxonomy table
# Expected schema:
#   Feature ID, Reference ID, Identity, AlignmentLength, Taxon
# Normalize the join key to "Sequence"
# ------------------------------------------------------------------------------
tax <- fread(tax_file, check.names = FALSE)

required_tax_cols <- c(
  "Feature ID",
  "Taxon",
  "Confidence",
  "Identity",
  "AlignmentLength",
  "TopHits"
)

missing_tax_cols <- setdiff(required_tax_cols, colnames(tax))

if (length(missing_tax_cols) > 0) {
  stop(sprintf(
    "VSEARCH taxonomy file is missing required columns: %s",
    paste(missing_tax_cols, collapse = ", ")
  ))
}


setnames(tax, "Feature ID", "Sequence")

tax[, Sequence := as.character(Sequence)]
tax[, Taxon := as.character(Taxon)]
tax[, Confidence := as.numeric(Confidence)]
tax[, Identity := as.numeric(Identity)]
tax[, AlignmentLength := as.numeric(AlignmentLength)]
tax[, TopHits := as.integer(TopHits)]

tax <- tax[, .(
  Sequence,
  Taxon,
  Confidence,
  Identity,
  AlignmentLength,
  TopHits
)]
# ------------------------------------------------------------------------------
# Read ASV abundance table
# Expected format:
#   first column = sample name
#   remaining columns = ASV sequences
# ------------------------------------------------------------------------------
asv <- fread(asv_file, check.names = FALSE)

if (ncol(asv) < 2) {
  stop("ASV table must contain at least one sample column and one ASV column")
}

setnames(asv, 1, "sample")

sample_names <- asv[["sample"]]
seq_cols <- setdiff(colnames(asv), "sample")

if (length(seq_cols) == 0) {
  stop("No ASV sequence columns found in ASV table")
}

# Convert sample x ASV table to ASV x sample table
abund_mat <- as.matrix(asv[, ..seq_cols])
rownames(abund_mat) <- sample_names

asv_long <- as.data.table(
  t(abund_mat),
  keep.rownames = "Sequence"
)

setnames(
  asv_long,
  old = setdiff(colnames(asv_long), "Sequence"),
  new = sample_names
)

# ------------------------------------------------------------------------------
# Merge VSEARCH taxonomy with abundance
# Keep all ASVs, including those without a VSEARCH hit
# ------------------------------------------------------------------------------
merged <- merge(
  asv_long,
  tax,
  by = "Sequence",
  all.x = TRUE,
  sort = FALSE
)

# Explicitly mark no-hit ASVs
merged[
  is.na(Taxon),
  Taxon := "Unclassified"
]

# ------------------------------------------------------------------------------
# Reorder columns
#
# ------------------------------------------------------------------------------
sample_cols <- setdiff(
  colnames(merged),
  c(
    "Sequence",
    "Taxon",
    "Confidence",
    "Identity",
    "AlignmentLength",
    "TopHits"
  )
)

setcolorder(
  merged,
  c(
    "Sequence",
    "Taxon",
    "Confidence",
    "Identity",
    "AlignmentLength",
    "TopHits",
    sample_cols
  )
)

fwrite(
  merged,
  out_file,
  sep = "\t",
  quote = FALSE,
  na = ""
)

message("Merged VSEARCH table written: ", out_file)
