#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

# Expect:
# 1) learn_nbases
# 2..n) filtered FASTQ files
if (length(args) < 3) {
  err_quit(
    paste(
      "Expected arguments:",
      "1) learn_nbases",
      "2) comma-separated binned quality scores, e.g. 3,10,17,22,27,35,40",
      "3..n) filtered FASTQ files",
      sep = "\n"
    )
  )
}

# --- parse arguments ---
learn_nbases <- suppressWarnings(as.numeric(args[1]))
if (is.na(learn_nbases) || learn_nbases <= 0) {
  err_quit("learn_nbases must be a positive number")
}

binnedQs <- suppressWarnings(as.numeric(strsplit(args[2], ",")[[1]]))

if (any(is.na(binnedQs)) || length(binnedQs) == 0) {
  err_quit("binned quality scores must be comma-separated numbers, e.g. 3,10,17,22,27,35,40")
}

filts <- args[-c(1, 2)]

# --- sanity checks ---
missing <- filts[!file.exists(filts)]
if (length(missing) > 0) {
  err_quit(
    paste(
      "Missing input FASTQ files:",
      paste(missing, collapse = "\n")
    )
  )
}

if (!exists("makeBinnedQualErrfun", where = asNamespace("dada2"), inherits = FALSE)) {
  err_quit("makeBinnedQualErrfun not found in dada2 namespace")
}

# --- build error function (PacBio CCS specific) ---
errfun <- dada2:::makeBinnedQualErrfun(binnedQs)

message("Learning error model")
message("Number of files: ", length(filts))
message("nbases: ", learn_nbases)

# --- learn errors ---
err <- suppressWarnings(
  learnErrors(
    fls = filts,
    nbases = learn_nbases,
    errorEstimationFunction = errfun,
    multithread = TRUE,
    randomize = TRUE,
    verbose = TRUE
  )
)

# --- outputs ---
saveRDS(err, "errorfun.rds")

pdf("plot_error_model.pdf", width = 12, height = 8, useDingbats = FALSE)
print(plotErrors(err))
dev.off()

message("Finished successfully")
quit(save = "no", status = 0)
