#!/usr/bin/env Rscript

###################################################
# This R script takes an input directory of .fastq.gz files
# and outputs a tsv file of the dada2 processed sequence
# table. It is intended for use with the QIIME2 plugin
# for DADA2.
#
# Ex: Rscript run_dada_single.R input_dir output.tsv track.tsv
#       primer_removed_dir filtered_dir front adapter 2 FALSE 0 0 2.0 2 20 Inf
#       pseudo pooled 1.0 0 1000000
####################################################

####################################################
#             DESCRIPTION OF ARGUMENTS             #
####################################################
# NOTE: All numeric arguments should be zero or positive.
# NOTE: All numeric arguments save maxEE are expected to be integers.
# NOTE: Currently the primer_removed_dir must already exist.
# NOTE: Currently the filterered_dir must already exist.
# NOTE: ALL ARGUMENTS ARE POSITIONAL!
#
### FILE SYSTEM ARGUMENTS ###
#
# 1) File path to directory with the .fastq.gz files to be processed.
#    Ex: path/to/dir/with/fastqgzs
#
# 2) File path to output tsv file. If already exists, will be overwritten.
#    Ex: path/to/output_file.tsv
#
# 3) File path to tracking tsv file. If already exists, will be overwritte.
#    Ex: path/to/tracking_stats.tsv
#
# 4) File path to directory in which to write the primer.removed .fastq.gz
#                 files. These files are intermediate for the full workflow.
#                 Currently they remain after the script finishes.
#                 Directory must already exist.
#    Ex: path/to/dir/with/fastqgzs/primerremoved
#
# 5) File path to directory in which to write the filtered .fastq.gz files.
#                 These files are intermediate for the full workflow.
#                 Currently they remain after the script finishes.
#                 Directory must already exist.
#    Ex: path/to/dir/with/fastqgzs/filtered
#
### PRIMER REMOVING ARGUMENTS ###
#
# 6) primerF - Primer front of Pacbio CCS sequences.
#    Ex: 'AGRGTTYGATYMTGGCTCAG'
#
# 7) primerR - Primer adapter of Pacbio CCS sequences.
#    Ex: 'RGYTACCTTGTTACGACTT'
#
# 8) maxMismatch - The number of mismatches to tolerate when matching
#                  reads to primer sequences.
#    Ex: 2
#
# 9) indels - Allow insertions or deletions of bases when matching adapters.
#    Ex: FALSE
#
### FILTERING ARGUMENTS ###
#
# 10) truncLen - The position at which to truncate reads. Reads shorter
#                than truncLen will be discarded.
#                Special values: 0 - no truncation or length filtering.
#     Ex: 150
#
# 11) trimLeft - The number of nucleotides to remove from the start of
#                each read. Should be less than truncLen for obvious reasons.
#     Ex: 0
#
# 12) maxEE - Reads with expected errors higher than maxEE are discarded.
#     Ex: 2.0
#
# 13) truncQ - Reads are truncated at the first instance of quality score truncQ.
#              If the read is then shorter than truncLen, it is discarded.
#     Ex: 2
#
# 14) minLen - Remove reads with length shorter than minLen. minLen is enforced
#              after trimming and truncation.
#              Default Inf - no maximum.
#     Ex: 20
#
# 15) maxLen - Remove reads with length greater than maxLen. maxLen is enforced
#              on the raw reads.
#              Default Inf - no maximum.
#    Ex: 300
#
### SENSITIVITY ARGUMENTS ###
#
# 16) poolMethod - The method used to pool (or not) samples during denoising.
#                  Valid options are:
#          independent: (Default) No pooling, samples are denoised indpendently.
#
#          pseudo: Samples are "pseudo-pooled" for denoising.
#    Ex: independent
#
#
### CHIMERA ARGUMENTS ###
#
# 17) chimeraMethod - The method used to remove chimeras. Valid options are:
#               none: No chimera removal is performed.
#               pooled: All reads are pooled prior to chimera detection.
#               consensus: Chimeras are detect in samples individually, and a
#                          consensus decision is made for each sequence variant.
#    Ex: consensus
#
# 18) minParentFold - The minimum abundance of potential "parents" of a sequence
#                     being tested as chimeric, expressed as a fold-change
#                     versus the abundance of the sequence being tested. Values
#                     should be greater than or equal to 1 (i.e. parents should
#                     be more abundant than the sequence being tested).
#    Ex: 1.0
#
### SPEED ARGUMENTS ###
#
# 19) nthreads - The number of threads to use.
#                 Special values: 0 - detect available cores and use all.
#    Ex: 1
#
# 20) nreads_learn - The minimum number of reads to learn the error model from.
#                 Special values: 0 - Use all input reads.
#    Ex: 1000000
#
#
#

cat(R.version$version.string, "\n")
errQuit <- function(mesg, status=1) { message("Error: ", mesg); q(status=status) }
args <- commandArgs(TRUE)

# Assign each of the arguments, in positional order, to an appropriately named R variable
inp.dir <- args[[1]]
out.path <- args[[2]]
out.track <- args[[3]]
primer.removed.dir <- args[[4]]
filtered.dir <- args[[5]]
primerF <- args[[6]]
primerR <- args[[7]]
maxMismatch <- as.numeric(args[[8]])
indels <- as.logical(args[[9]])
truncLen <- as.integer(args[[10]])
trimLeft <- as.integer(args[[11]])
maxEE <- as.numeric(args[[12]])
truncQ <- as.integer(args[[13]])
minLen <- as.numeric(args[[14]])
maxLen <- as.numeric(args[[15]]) # Allows Inf
poolMethod <- args[[16]]
chimeraMethod <- args[[17]]
minParentFold <- as.numeric(args[[18]])
nthreads <- as.integer(args[[19]])
nreads.learn <- as.integer(args[[20]])

### VALIDATE ARGUMENTS ###

# Input directory is expected to contain .fastq.gz file(s)
# that have not yet been filtered and globally trimmed
# to the same length.
if(!dir.exists(inp.dir)) {
  errQuit("Input directory does not exist.")
} else {
  unfilts <- list.files(inp.dir, pattern=".fastq.gz$", full.names=TRUE)
  if(length(unfilts) == 0) {
    errQuit("No input files with the expected filename format found.")
  }
}

# Output files are to be filenames (not directories) and are to be
# removed and replaced if already present.
for(fn in c(out.path, out.track)) {
  if(dir.exists(fn)) {
    errQuit("Output filename ", fn, " is a directory.")
  } else if(file.exists(fn)) {
    invisible(file.remove(fn))
  }
}

# Convert nthreads to the logical/numeric expected by dada2
if(nthreads < 0) {
  errQuit("nthreads must be non-negative.")
} else if(nthreads == 0) {
  multithread <- TRUE # detect and use all
} else if(nthreads == 1) {
  multithread <- FALSE
} else {
  multithread <- nthreads
}

### LOAD LIBRARIES ###
suppressWarnings(library(methods))
suppressWarnings(library(dada2))
cat("DADA2:", as.character(packageVersion("dada2")), "/",
    "Rcpp:", as.character(packageVersion("Rcpp")), "/",
    "RcppParallel:", as.character(packageVersion("RcppParallel")), "\n")

### Remove Primers ###
cat("1) Removing Primers\n")
nop <- file.path(primer.removed.dir, basename(unfilts))
if(primerF == 'none' && primerR == 'none'){
  nop <- unfilts
} else {
  prim <- suppressWarnings(removePrimers(unfilts, nop, primerF, dada2:::rc(primerR),
                                         max.mismatch = maxMismatch, allow.indels = indels,
                                         orient = TRUE, verbose = TRUE))
  cat(ifelse(file.exists(nop), ".", "x"), sep="")
  nop <- list.files(primer.removed.dir, pattern=".fastq.gz$", full.names=TRUE)
  cat("\n")
  if(length(nop) == 0) { # All reads were filtered out
    errQuit("No reads passed the Removing Primers step  (Did you select the right primers?)", status=2)
  }
}

### TRIM AND FILTER ###
cat("2) Filtering\n")
filts <- file.path(filtered.dir, basename(nop))
out <- suppressWarnings(filterAndTrim(nop, filts, truncLen = truncLen, trimLeft = trimLeft,
                                      maxEE = maxEE, truncQ = truncQ, rm.phix = FALSE,
                                      multithread = multithread, maxLen = maxLen, minLen = minLen, minQ = 0))
cat(ifelse(file.exists(filts), ".", "x"), sep="")
filts <- list.files(filtered.dir, pattern=".fastq.gz$", full.names=TRUE)
cat("\n")
if(length(filts) == 0) { # All reads were filtered out
  errQuit("No reads passed the filter (was truncLen longer than the read length?)", status=2)
}

### LEARN ERROR RATES ###
# Dereplicate enough samples to get nreads.learn total reads
cat("3) Learning Error Rates\n")
err <- suppressWarnings(learnErrors(filts, nreads=nreads.learn,
                                    errorEstimationFunction=dada2:::PacBioErrfun,
                                    multithread=multithread, BAND_SIZE=32))

### PROCESS ALL SAMPLES ###
# Loop over rest in streaming fashion with learned error rates
dds <- vector("list", length(filts))
cat("4) Denoise samples ")
cat("\n")
for(j in seq(length(filts))) {
  drp <- derepFastq(filts[[j]])
  dds[[j]] <- dada(drp, err=err, multithread=multithread,
                   BAND_SIZE=32, verbose=FALSE)
  cat(".")
}
cat("\n")
if(poolMethod == "pseudo") {
  cat("  Pseudo-pool step ")
  ### TEMPORARY, to be removed once 1.12 makes its way to Q2
  ### Needed for now to manage pseudo-pooling memory, as 1.10 didn't do this appropriately.
  ### pseudo_priors code copied from dada2.R
  st <- makeSequenceTable(dds)
  pseudo_priors <- colnames(st)[colSums(st>0) >= 2 | colSums(st) >= Inf]
  rm(st)
  ### \pseudo_priors code copied from dada2.R
  ### code copied from previous loop through samples in this script
  for(j in seq(length(filts))) {
    drp <- derepFastq(filts[[j]])
    dds[[j]] <- dada(drp, err=err, multithread=multithread,
                     priors = pseudo_priors,
                     BAND_SIZE=32, verbose=FALSE)
    cat(".")
  }
  cat("\n")
  ### \code copied from previous loop through samples in this script
}

# Make sequence table
seqtab <- makeSequenceTable(dds)

### Remove chimeras
cat("5) Remove chimeras (method = ", chimeraMethod, ")\n", sep="")
if(chimeraMethod %in% c("pooled", "consensus")) {
  seqtab.nochim <- removeBimeraDenovo(seqtab, method=chimeraMethod, minFoldParentOverAbundance=minParentFold, multithread=multithread)
} else { # No chimera removal, copy seqtab to seqtab.nochim
  seqtab.nochim <- seqtab
}

### REPORT READ FRACTIONS THROUGH PIPELINE ###
cat("6) Report read numbers through the pipeline\n")
# Handle edge cases: Samples lost in filtering; One sample
if (primerF == 'none' && primerR == 'none'){
  track <- cbind(out[ ,1:2,drop=FALSE], matrix(0, nrow=nrow(out), ncol=2))
  colnames(track) <- c("input", "filtered", "denoised", "non-chimeric")
  passed.filtering <- track[,"filtered"] > 0
  track[passed.filtering,"denoised"] <- rowSums(seqtab)
  track[passed.filtering,"non-chimeric"] <- rowSums(seqtab.nochim)
  write.table(track, out.track, sep="\t", row.names=TRUE, col.names=NA,
              quote=FALSE)
} else {
  track <- cbind(prim,out[ ,2], matrix(0, nrow=nrow(out), ncol=2))
  colnames(track) <- c("input", "primer-removed","filtered", "denoised", "non-chimeric")
  passed.filtering <- track[,"filtered"] > 0
  track[passed.filtering,"denoised"] <- rowSums(seqtab)
  track[passed.filtering,"non-chimeric"] <- rowSums(seqtab.nochim)
  write.table(track, out.track, sep="\t", row.names=TRUE, col.names=NA,
              quote=FALSE)
}

### WRITE OUTPUT AND QUIT ###
# Formatting as tsv plain-text sequence table table
cat("7) Write output\n")
seqtab.nochim <- t(seqtab.nochim) # QIIME has OTUs as rows
col.names <- basename(filts)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab.nochim, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
saveRDS(seqtab.nochim, gsub("tsv", "rds", out.path)) ### TESTING
saveRDS(seqtab.nochim, "./seqtab_nochim.rds") # Added by Khi Pin

q(status=0)
