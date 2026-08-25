#!/usr/bin/env Rscript

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 8) {
  err_quit(
    paste(
      "Expected at least 8 positional arguments:",
      "1) input_seqtab_nochim_rds",
      "2) input_seqtab_filtered_rds",
      "3) metadata_tsv",
      "4) output_dada2_tracking_tsv",
      "5) output_frequency_tsv",
      "6) output_rarefaction_depth_txt",
      "7..n) filter_stats files, then literal --DENOISE_STATS--, then denoise_stats files",
      sep = "\n"
    )
  )
}

marker_idx <- match("--DENOISE_STATS--", args)

if (is.na(marker_idx)) {
  err_quit("Missing required separator token: --DENOISE_STATS--")
}

if (marker_idx <= 7) {
  err_quit("No filter_stats files provided before --DENOISE_STATS--")
}

if (marker_idx == length(args)) {
  err_quit("No denoise_stats files provided after --DENOISE_STATS--")
}


# -------------------------------------------------------------------------
# Arguments
# -------------------------------------------------------------------------

input_seqtab_nochim_rds <- args[1]
input_seqtab_filtered_rds <- args[2]
metadata_tsv <- args[3]

output_tracking_tsv <- args[4]
output_frequency_tsv <- args[5]
output_rarefaction_txt <- args[6]

filter_stats_files <- args[7:(marker_idx - 1)]
denoise_stats_files <- args[(marker_idx + 1):length(args)]


# -------------------------------------------------------------------------
# Validate input files
# -------------------------------------------------------------------------

for (f in c(
  input_seqtab_nochim_rds,
  input_seqtab_filtered_rds,
  metadata_tsv,
  filter_stats_files,
  denoise_stats_files
)) {
  if (!file.exists(f)) {
    err_quit(sprintf("Required file does not exist: %s", f))
  }
}


# -------------------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------------------

read_tsv_safe <- function(path) {
  read.table(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


bind_named_tables <- function(files, required_cols) {
  tabs <- lapply(files, read_tsv_safe)

  for (i in seq_along(tabs)) {
    missing_cols <- setdiff(required_cols, colnames(tabs[[i]]))

    if (length(missing_cols) > 0) {
      err_quit(
        sprintf(
          "File %s is missing required columns: %s",
          files[i],
          paste(missing_cols, collapse = ", ")
        )
      )
    }
  }

  out <- do.call(rbind, tabs)
  rownames(out) <- NULL

  out
}


# -------------------------------------------------------------------------
# Read filtering and denoising statistics
# -------------------------------------------------------------------------

filter_stats <- bind_named_tables(
  filter_stats_files,
  c("sample", "input_reads", "filtered_reads")
)

denoise_stats <- bind_named_tables(
  denoise_stats_files,
  c("sample", "denoised_reads")
)


# -------------------------------------------------------------------------
# Read DADA2 sequence tables
# -------------------------------------------------------------------------

seqtab_nochim <- readRDS(input_seqtab_nochim_rds)
seqtab_filtered <- readRDS(input_seqtab_filtered_rds)

if (!is.matrix(seqtab_nochim)) {
  err_quit("seqtab_nochim is not a matrix")
}

if (!is.matrix(seqtab_filtered)) {
  err_quit("seqtab_filtered is not a matrix")
}


# -------------------------------------------------------------------------
# DADA2 read tracking
# -------------------------------------------------------------------------

nochim_df <- data.frame(
  sample = rownames(seqtab_nochim),
  non_chimeric_reads = as.numeric(rowSums(seqtab_nochim)),
  stringsAsFactors = FALSE
)

tracking <- merge(
  filter_stats,
  denoise_stats,
  by = "sample",
  all = TRUE,
  sort = FALSE
)

tracking <- merge(
  tracking,
  nochim_df,
  by = "sample",
  all = TRUE,
  sort = FALSE
)

tracking <- tracking[
  ,
  c(
    "sample",
    "input_reads",
    "filtered_reads",
    "denoised_reads",
    "non_chimeric_reads"
  )
]

sample_order <- unique(
  c(
    filter_stats$sample,
    denoise_stats$sample,
    nochim_df$sample
  )
)

tracking <- tracking[
  match(sample_order, tracking$sample),
]

tracking <- tracking[!is.na(tracking$sample), ]

write.table(
  tracking,
  file = output_tracking_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# -------------------------------------------------------------------------
# Final read abundance and observed ASV richness
# -------------------------------------------------------------------------

frequency_df <- data.frame(
  sample = rownames(seqtab_filtered),
  reads_after_asv_filtering = as.numeric(rowSums(seqtab_filtered)),
  number_of_asvs = as.numeric(rowSums(seqtab_filtered > 0)),
  stringsAsFactors = FALSE
)


# -------------------------------------------------------------------------
# Add sample metadata
# -------------------------------------------------------------------------

metadata <- read_tsv_safe(metadata_tsv)

if (!("sample_name" %in% colnames(metadata))) {
  err_quit("Metadata must contain column 'sample_name'")
}

frequency_df <- merge(
  frequency_df,
  metadata,
  by.x = "sample",
  by.y = "sample_name",
  all.x = TRUE,
  sort = FALSE
)

frequency_df <- frequency_df[
  match(rownames(seqtab_filtered), frequency_df$sample),
]

write.table(
  frequency_df,
  file = output_frequency_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)


# -------------------------------------------------------------------------
# Suggested rarefaction depth
# -------------------------------------------------------------------------

sorted_depths <- sort(
  frequency_df$reads_after_asv_filtering,
  decreasing = TRUE
)

number <- length(sorted_depths)

if (number == 0) {
  err_quit("No samples available in final filtered sequence table")
}

if (number <= 2) {

  rarefaction_d <- sorted_depths[number]

} else {

  # Select a depth reached by approximately 80% of samples.
  result <- as.integer(number * 0.8)

  if (result < 1) {
    result <- 1L
  }

  rarefaction_d <- sorted_depths[result]
}

writeLines(
  as.character(as.integer(floor(rarefaction_d))),
  output_rarefaction_txt
)

quit(save = "no", status = 0)
