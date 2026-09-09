#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(
    paste(
      "Usage:",
      "merge_vsearch_taxonomy_with_table.R",
      "<vsearch_lca.tsv>",
      "<vsearch_hits.tsv>",
      "<asv_table.tsv>",
      "<out.tsv>"
    )
  )
}

tax_file <- args[1]
hits_file <- args[2]
asv_file <- args[3]
out_file <- args[4]

for (f in c(tax_file, hits_file, asv_file)) {
  if (!file.exists(f)) {
    stop(sprintf("Input file not found: %s", f))
  }
}

# ------------------------------------------------------------------------------
# Read native VSEARCH --lcaout output
#
# Expected format:
#
#   <Sequence>    <Taxon>
#
# No header.
# Unassigned sequences may have an empty taxonomy field.
# ------------------------------------------------------------------------------

tax <- fread(
  tax_file,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  check.names = FALSE
)

if (ncol(tax) < 1) {
  stop("VSEARCH LCA file contains no columns")
}

if (ncol(tax) == 1) {
  tax[, V2 := NA_character_]
}

tax <- tax[, .(
  Sequence = as.character(V1),
  Taxon = as.character(V2)
)]

tax[
  is.na(Taxon) | trimws(Taxon) == "",
  Taxon := NA_character_
]

# Normalize any trailing semicolon so taxonomy comparisons are consistent.
tax[
  !is.na(Taxon),
  Taxon := sub(";$", "", Taxon)
]

if (anyDuplicated(tax$Sequence)) {
  stop("VSEARCH LCA file contains duplicate query sequences")
}

# ------------------------------------------------------------------------------
# Read VSEARCH --userout hits
#
# Expected userfields:
#
#   query+target+id+qcov+alnlen+ql+tl
#
# Therefore:
#
#   V1 = query
#   V2 = target
#   V3 = identity
#   V4 = query coverage
#   V5 = alignment length
#   V6 = query length
#   V7 = target length
# ------------------------------------------------------------------------------

hits <- fread(
  hits_file,
  header = FALSE,
  sep = "\t",
  fill = TRUE,
  quote = "",
  check.names = FALSE
)

if (nrow(hits) > 0) {
  if (ncol(hits) < 7) {
    stop(
      sprintf(
        "VSEARCH hits file has %d columns; expected at least 7",
        ncol(hits)
      )
    )
  }

  hits <- hits[, .(
    Sequence = as.character(V1),
    Target = as.character(V2),
    Identity = as.numeric(V3),
    QueryCoverage = as.numeric(V4),
    AlignmentLength = as.numeric(V5),
    QueryLength = as.numeric(V6),
    TargetLength = as.numeric(V7)
  )]

  # ---------------------------------------------------------------------------
  # Extract taxonomy from SINTAX-style target headers
  #
  # Example:
  #
  #   KY103007;tax=d:Fungi,p:Basidiomycota,...,s:dermatis;
  #
  # becomes:
  #
  #   d:Fungi,p:Basidiomycota,...,s:dermatis
  # ---------------------------------------------------------------------------

  hits[, HitTaxon := sub("^.*;tax=", "", Target)]
  hits[, HitTaxon := sub(";$", "", HitTaxon)]

  # Targets without a tax= field are not useful for support calculation.
  hits[!grepl(";tax=", Target, fixed = TRUE), HitTaxon := NA_character_]

} else {

  hits <- data.table(
    Sequence = character(),
    Target = character(),
    Identity = numeric(),
    QueryCoverage = numeric(),
    AlignmentLength = numeric(),
    QueryLength = numeric(),
    TargetLength = numeric(),
    HitTaxon = character()
  )
}

# ------------------------------------------------------------------------------
# Calculate support for the LCA assignment
#
# Confidence is defined as:
#
#   percentage of accepted VSEARCH hits supporting the complete lineage
#   reported by --lcaout.
#
# Example:
#
# LCA:
#   d:Fungi,p:Basidiomycota,...,g:Cutaneotrichosporon
#
# Hits:
#   ...,g:Cutaneotrichosporon,s:dermatis
#   ...,g:Cutaneotrichosporon,s:faecale
#   ...,g:Cutaneotrichosporon,s:mucoides
#
# All support the assigned genus, therefore Confidence = 100.
#
# This is NOT equivalent to DADA2 bootstrap confidence.
# ------------------------------------------------------------------------------

calculate_support <- function(sequence, assigned_taxon) {

  if (
    is.na(assigned_taxon) ||
    assigned_taxon == ""
  ) {
    return(NA_real_)
  }

  query_hits <- hits[
    Sequence == sequence &
      !is.na(HitTaxon)
  ]

  if (nrow(query_hits) == 0) {
    return(NA_real_)
  }

  # A hit supports the assignment when its taxonomy begins with the full
  # LCA lineage and either ends there or continues to a lower rank.
  #
  # Example:
  #
  # assigned:
  #   d:Fungi,p:Ascomycota
  #
  # supported:
  #   d:Fungi,p:Ascomycota,c:Eurotiomycetes
  #
  # not supported:
  #   d:Fungi,p:Basidiomycota
  prefix <- paste0(assigned_taxon, ",")

  supported <- (
    query_hits$HitTaxon == assigned_taxon |
      startsWith(query_hits$HitTaxon, prefix)
  )

  100 * sum(supported) / length(supported)
}

tax[, Confidence := mapply(
  calculate_support,
  Sequence,
  Taxon
)]

# Keep confidence consistent with the integer 0-100 convention used elsewhere.
tax[
  !is.na(Confidence),
  Confidence := round(Confidence)
]

# ------------------------------------------------------------------------------
# Read ASV abundance table
#
# Expected format:
#
#   first column      = sample name
#   remaining columns = ASV sequences
# ------------------------------------------------------------------------------

asv <- fread(
  asv_file,
  check.names = FALSE
)

if (ncol(asv) < 2) {
  stop(
    "ASV table must contain at least one sample column and one ASV column"
  )
}

setnames(asv, 1, "sample")

sample_names <- as.character(asv[["sample"]])
seq_cols <- setdiff(colnames(asv), "sample")

if (length(seq_cols) == 0) {
  stop("No ASV sequence columns found in ASV table")
}

if (anyDuplicated(sample_names)) {
  stop("ASV table contains duplicate sample names")
}

# ------------------------------------------------------------------------------
# Convert sample x ASV table to ASV x sample table
# ------------------------------------------------------------------------------

abund_mat <- as.matrix(
  asv[, ..seq_cols]
)

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
# Merge taxonomy with abundance
#
# Keep every ASV from the abundance table, even if VSEARCH did not assign it.
# ------------------------------------------------------------------------------

merged <- merge(
  asv_long,
  tax,
  by = "Sequence",
  all.x = TRUE,
  sort = FALSE
)

# ------------------------------------------------------------------------------
# Final column order
#
# Sequence, Taxon, Confidence, <sample columns...>
#
# The generic taxonomy_add_md5 process will prepend the final "id" column.
# ------------------------------------------------------------------------------

sample_cols <- setdiff(
  colnames(merged),
  c("Sequence", "Taxon", "Confidence")
)

setcolorder(
  merged,
  c(
    "Sequence",
    "Taxon",
    "Confidence",
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

message("Merged VSEARCH taxonomy table written: ", out_file)
