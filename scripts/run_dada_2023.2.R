#!/usr/bin/env Rscript

###################################################
# This R script takes an input directory of .fastq.gz files
# and outputs a tsv file of the dada2 processed sequence
# table. It is intended for use with the QIIME2 plugin
# for DADA2.
#
# Ex: Rscript run_dada_single.R input_dir output.tsv track.tsv filtered_dir 200 0 2.0 2 Inf pooled 1.0 0 1000000 NULL 32
####################################################

####################################################
#             DESCRIPTION OF ARGUMENTS             #
####################################################
# NOTE: All numeric arguments should be zero or positive.
# NOTE: All numeric arguments save maxEE are expected to be integers.
# NOTE: Currently the filterered_dir must already exist.
# NOTE: ALL ARGUMENTS ARE POSITIONAL!
#
### FILE SYSTEM ARGUMENTS ###
#
# 1) input_directory - File path to directory with the .fastq.gz files to be processed.
#    Ex: path/to/dir/with/fastqgzs
#
# 2) output_path - File path to output tsv file. If already exists, will be overwritten.
#    Ex: path/to/output_file.tsv
#
# 3) output_track - File path to tracking tsv file. If already exists, will be overwritten.
#    Ex: path/to/tracking_stats.tsv
#
# 4) removed_primer_directory - File path to directory in which to write the primer.removed .fastq.gz
#                 files. These files are intermediate for the full workflow.
#                 Currently they remain after the script finishes.
#                 Directory must already exist.
#    Ex: path/to/dir/with/fastqgzs/primerremoved
#
# 5) filtered_directory - File path to directory in which to write the filtered .fastq.gz files.
#                 These files are intermediate for the full workflow.
#                 Currently they remain after the script finishes.
#                 Directory must already exist.
#    Ex: path/to/dir/with/fastqgzs/filtered
#
### PRIMER REMOVING ARGUMENTS ###
#
# 6) forward_primer - Primer front of Pacbio CCS sequences.
#    Ex: 'AGRGTTYGATYMTGGCTCAG'
#
# 7) reverse_primer - Primer adapter of Pacbio CCS sequences.
#    Ex: 'RGYTACCTTGTTACGACTT'
#
# 8) max_mismatch - The number of mismatches to tolerate when matching
#                  reads to primer sequences.
#    Ex: 2
#
# 9) indels - Allow insertions or deletions of bases when matching adapters.
#    Ex: FALSE
#
### FILTERING ARGUMENTS ###
#
# 10) truncation_length - The position at which to truncate reads. Reads shorter
#               than truncation_length  will be discarded.
#               Special values: 0 - no truncation or length filtering.
#    Ex: 150
#
# 11) trim_left - The number of nucleotides to remove from the start of
#               each read. Should be less than truncation_length  for obvious reasons.
#    Ex: 0
#
# 12) max_expected_errors - Reads with expected errors higher than maxEE are discarded.
#    Ex: 2.0
#
# 13) truncation_quality_score - Reads are truncated at the first instance of quality score truncation_quality_score.
#                If the read is then shorter than truncation_length , it is discarded.
#    Ex: 2
#
# 14) min_length - Remove reads with length shorter than min_length. min_length is enforced
#              after trimming and truncation.
#              Default Inf - no maximum.
#
# 15) max_length - Remove reads with length greater than max_length. max_length is enforced on the raw reads.
#             Default Inf - no maximum.
#    Ex: 300
#
### SENSITIVITY ARGUMENTS ###
#
# 16) pooling_method- The method used to pool (or not) samples during denoising.
#             Valid options are:
#               independent: (Default) No pooling, samples are denoised indpendently.
#               pseudo: Samples are "pseudo-pooled" for denoising.
#    Ex: independent
#
#
### CHIMERA ARGUMENTS ###
#
# 17) chimera_method - The method used to remove chimeras. Valid options are:
#               none: No chimera removal is performed.
#               pooled: All reads are pooled prior to chimera detection.
#               consensus: Chimeras are detected in samples individually, and a consensus decision
#                           is made for each sequence variant.
#    Ex: consensus
#
# 18) min_parental_fold - The minimum abundance of potential "parents" of a sequence being
#               tested as chimeric, expressed as a fold-change versus the abundance of the sequence being
#               tested. Values should be greater than or equal to 1 (i.e. parents should be more
#               abundant than the sequence being tested).
#    Ex: 1.0
#
### SPEED ARGUMENTS ###
#
# 19) num_threads - The number of threads to use.
#                 Special values: 0 - detect available cores and use all.
#    Ex: 1
#
# 20) learn_min_reads - The minimum number of reads to learn the error model from.
#                 Special values: 0 - Use all input reads.
#    Ex: 1000000
#
### GLOBAL OPTION ARGUMENTS ###
#
# 21) homopolymer_gap_penalty - The cost of gaps in homopolymer regions (>=3 repeated bases).
#                               Default is NULL, which causes homopolymer gaps
#                               to be treated as normal gaps.
#    Ex: -1
#
# 22) band_size - When set, banded Needleman-Wunsch alignments are performed.
#                 The default value of band_size is 16. Setting BAND_SIZE to a negative
#                 number turns off banding (i.e. full Needleman-Wunsch).
#    Ex: 32
#

# error handling -----------------
options(error = function() {
  sink(stderr())
  on.exit(sink(NULL))
  traceback(3)
  if (!interactive()) {
    q(status = 1)
  }
})

library("optparse")

cat(R.version$version.string, "\n")
errQuit <- function(mesg, status=1) { message("Error: ", mesg); q(status=status) }
getN <- function(x) sum(getUniques(x)) #Function added from paired read processing

option_list = list(
  make_option(c("--input_directory"), action="store", default='NULL', type='character',
              help="File path to directory with the .fastq.gz files to be processed"),
  make_option(c("--input_directory_reverse"), action="store", default='NULL', type='character',
              help="File path to reverse directory with the .fastq.gz files to be processed. Only used in paired-end processing"),
  make_option(c("--output_path"), action="store", default='NULL', type='character',
              help="File path to output tsv file. If already exists, will be overwritten"),
  make_option(c("--output_track"), action="store", default='NULL', type='character',
              help="File path to tracking tsv file. If already exists, will be overwritten"),
  make_option(c("--removed_primer_directory"), action="store", default='NULL', type='character',
              help="File path to directory in which to write the primer.removed .fastq.gz files"),
  make_option(c("--filtered_directory"), action="store", default='NULL', type='character',
              help="File path to directory in which to write the filtered .fastq.gz files. These files are intermediate"),
  make_option(c("--filtered_directory_reverse"), action="store", default='NULL', type='character',
              help="File path to directory in which to write the reverse filtered .fastq.gz files. These files are intermediate. Only used in paired-end processing"),
  make_option(c("--forward_primer"), action="store", default='NULL', type='character',
              help="Primer front of Pacbio CCS sequences"),
  make_option(c("--reverse_primer"), action="store", default='NULL', type='character',
              help="Primer adapter of Pacbio CCS sequences"),
  make_option(c("--max_mismatch"), action="store", default='NULL', type='character',
              help="The number of mismatches to tolerate when matching reads to primer sequences."),
  make_option(c("--indels"), action="store", default='NULL', type='character',
              help="Allow insertions or deletions of bases when matching adapters"),
  make_option(c("--truncation_length"), action="store", default='NULL', type='character',
              help="The position at which to truncate reads. Reads shorter then truncation_length will be discarded."),
  make_option(c("--truncation_length_reverse"), action="store", default='NULL', type='character',
              help="The position at which to truncate reverse reads. Reads shorter then truncation_length will be discarded. Only used in paired-end processing"),
  make_option(c("--trim_left"), action="store", default='NULL', type='character',
              help="The number of nucleotides to remove from the start of each read. Should be less than truncation_length for obvious reasons"),
  make_option(c("--trim_left_reverse"), action="store", default='NULL', type='character',
              help="The number of nucleotides to remove from the start of each reverse read. Should be less than truncation_length for obvious reasons. Only used in paired-end processing"),
  make_option(c("--max_expected_errors"), action="store", default='NULL', type='character',
              help="Reads with expected errors higher than max_expected_errors are discarded"),
  make_option(c("--max_expected_errors_reverse"), action="store", default='NULL', type='character',
              help="Reverse reads with expected errors higher than max_expected_errors are discarded. Only used in paired-end processing"),
  make_option(c("--truncation_quality_score"), action="store", default='NULL', type='character',
              help="Reads are truncated at the first instance of quality score truncation_quality_score.If the read is then shorter than truncation_length, it is discarded"),
  make_option(c("--min_length"), action="store", default='NULL', type='character',
              help="Remove reads with length shorter than min_length. min_length is enforced after trimming and truncation."),
  make_option(c("--max_length"), action="store", default='NULL', type='character',
              help="Remove reads with length greater than max_length. max_length is enforced on the raw reads."),
  make_option(c("--min_overlap"), action="store", default='NULL', type='character',
              help="Minimum overlap allowed between reads"),
  make_option(c("--pooling_method"), action="store", default='NULL', type='character',
              help="The method used to pool (or not) samples during denoising (independent/pseudo)"),
  make_option(c("--chimera_method"), action="store", default='NULL', type='character',
              help="The method used to remove chimeras (none/pooled/consensus)"),
  make_option(c("--min_parental_fold"), action="store", default='NULL', type='character',
              help="The minimum abundance of potential parents of a sequence being tested as chimeric, expressed as a fold-change versus the abundance of the sequence being tested. Values should be greater than or equal to 1"),
  make_option(c("--allow_one_off"), action="store", default='NULL', type='character',
              help="Bimeras that are one-off (one mismatch/indel away from an exact bimera) are also identified if allow_one_off is TRUE."),
  make_option(c("--num_threads"), action="store", default='NULL', type='character',
              help="The number of threads to use"),
  make_option(c("--learn_min_reads"), action="store", default='NULL', type='character',
              help="The minimum number of reads to learn the error model from"),
  make_option(c("--homopolymer_gap_penalty"), action="store", default='NULL', type='character',
              help="The cost of gaps in homopolymer regions (>=3 repeated bases).Default is NULL, which causes homopolymer gaps to be treated as normal gaps."),
  make_option(c("--band_size"), action="store", default='NULL', type='character',
              help="When set, banded Needleman-Wunsch alignments are performed.")
)
opt = parse_args(OptionParser(option_list=option_list))

# Added by Khi Pin 2023-4-21 to save all outputs from scripts
con <- file("dada2.log")
sink(con, append=TRUE)
sink(con, append=TRUE, type="message")


# Assign each of the arguments, in positional order, to an appropriately named R variable
inp.dir <- opt$input_directory
inp.dirR <- opt$input_directory_reverse #added from paired arguments
out.path <- opt$output_path
out.track <- opt$output_track
primer.removed.dir <- opt$removed_primer_directory #added from CCS arguments
filtered.dir <- opt$filtered_directory
filtered.dirR<- opt$filtered_directory_reverse #added from paired arguments
primerF <- opt$forward_primer #added from CCS arguments
primerR <- opt$reverse_primer #added from CCS arguments
maxMismatch <- if(opt$max_mismatch=='NULL') NULL else as.numeric(opt$max_mismatch) #added from CCS arguments
indels <- if(opt$indels=='NULL') NULL else as.logical(opt$indels)  #added from CCS arguments
truncLen <- if(opt$truncation_length=='NULL') NULL else as.integer(opt$truncation_length)
truncLenR <- if(opt$truncation_length_reverse=='NULL') NULL else as.integer(opt$truncation_length_reverse) #added from paired arguments
trimLeft <- if(opt$trim_left=='NULL') NULL else as.integer(opt$trim_left)
trimLeftR<-if(opt$trim_left_reverse=='NULL') NULL else as.integer(opt$trim_left_reverse) #added from paired arguments
maxEE <- if(opt$max_expected_errors=='NULL') NULL else as.numeric(opt$max_expected_errors)
maxEER <- if(opt$max_expected_errors_reverse=='NULL') NULL else as.numeric(opt$max_expected_errors_reverse) #added from paired arguments
truncQ <- if(opt$truncation_quality_score=='NULL') NULL else as.integer(opt$truncation_quality_score)
minLen <- if(opt$min_length=='NULL') NULL else as.numeric(opt$min_length) #added from CCS arguments
maxLen <- if(opt$max_length=='NULL') NULL else as.numeric(opt$max_length) # Allows Inf
minOverlap <- if(opt$min_overlap=='NULL') NULL else as.integer(opt$min_overlap) #added from paired arguments
poolMethod <- opt$pooling_method
chimeraMethod <- opt$chimera_method
minParentFold <- if(opt$min_parental_fold=='NULL') NULL else as.numeric(opt$min_parental_fold)
allowOneOff <-if(opt$allow_one_off=='NULL') NULL else as.logical(opt$allow_one_off)
nthreads <- if(opt$num_threads=='NULL') NULL else as.integer(opt$num_threads)
nreads.learn <- if(opt$learn_min_reads=='NULL') NULL else as.integer(opt$learn_min_reads)
# The following args are not directly exposed to end users in q2-dada2,
# but rather indirectly, via the methods `denoise-single` and `denoise-pyro`.
if (opt$homopolymer_gap_penalty=='NULL'){
  HOMOPOLYMER_GAP_PENALTY<-NULL
}else{
  HOMOPOLYMER_GAP_PENALTY<-as.integer(opt$homopolymer_gap_penalty)
  if(HOMOPOLYMER_GAP_PENALTY>0){
    HOMOPOLYMER_GAP_PENALTY<-HOMOPOLYMER_GAP_PENALTY*(-1) #negative numbers cannot be passed using optparse, so we convert it here
  }
}
BAND_SIZE <- if(opt$band_size=='NULL') NULL else as.integer(opt$band_size)

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
  if(inp.dirR!='NULL'){
    unfiltsR <- list.files(inp.dirR, pattern=".fastq.gz$", full.names=TRUE)
    if(length(unfiltsR) == 0) {
      errQuit("No input reverse files with the expected filename format found.")
    }
    if(length(unfilts) != length(unfiltsR)) {
      errQuit("Different numbers of forward and reverse .fastq.gz files.")
    }

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
if(primer.removed.dir!='NULL'){ #for CCS read analysis
  cat("1) Removing Primers\n")
  nop <- file.path(primer.removed.dir, basename(unfilts))

  if(primerF == 'none' && primerR == 'none'){
    nop <- unfilts
  } else {
    prim <- suppressWarnings(removePrimers(unfilts, nop, primer, dada2:::rc(primerR),
                                           max.mismatch = maxMismatch, allow.indels = indels,
                                           orient = TRUE, verbose = TRUE))
    cat(ifelse(file.exists(nop), ".", "x"), sep="")
    nop <- list.files(primer.removed.dir, pattern=".fastq.gz$", full.names=TRUE)
    cat("\n")
    if(length(nop) == 0) { # All reads were filtered out
      errQuit("No reads passed the Removing Primers step  (Did you select the right primers?)", status=2)
    }
  }
}

### TRIM AND FILTER ###
cat("2) Filtering ")
if(primer.removed.dir!='NULL'){ #for CCS read analysis
  filts <- file.path(filtered.dir, basename(nop))
  out <- suppressWarnings(filterAndTrim(nop, filts, truncLen = truncLen, trimLeft = trimLeft,
                                        maxEE = maxEE, truncQ = truncQ, rm.phix = FALSE,
                                        multithread = multithread, maxLen = maxLen, minLen = minLen, minQ = 0))
}else{
  filts <- file.path(filtered.dir, basename(unfilts))
  if(inp.dirR!='NULL'){#for paired read analysis
    filtsR <- file.path(filtered.dirR, basename(unfiltsR))
    out <- suppressWarnings(filterAndTrim(unfilts, filts, unfiltsR, filtsR,
                                          truncLen=c(truncLen, truncLenR), trimLeft=c(trimLeft, trimLeftR),
                                          maxEE=c(maxEE, maxEER), truncQ=truncQ, rm.phix=TRUE,
                                          multithread=multithread))
  }else{#for sinlge/pyro read analysis
    out <- suppressWarnings(filterAndTrim(unfilts, filts, truncLen=truncLen, trimLeft=trimLeft,
                                          maxEE=maxEE, truncQ=truncQ, rm.phix=TRUE,
                                          multithread=multithread, maxLen=maxLen))
  }
}

cat(ifelse(file.exists(filts), ".", "x"), sep="")
filts <- list.files(filtered.dir, pattern=".fastq.gz$", full.names=TRUE)

if(inp.dirR!='NULL'){
  filtsR <- list.files(filtered.dirR, pattern=".fastq.gz$", full.names=TRUE)
}
cat("\n")
if(length(filts) == 0) { # All reads were filtered out
  errQuit("No reads passed the filter (was truncLen longer than the read length?)", status=2)
}

### LEARN ERROR RATES ###
# Dereplicate enough samples to get nreads.learn total reads
cat("3) Learning Error Rates\n")
if(primer.removed.dir!='NULL'){#for CCS read analysis
  err <- suppressWarnings(learnErrors(filts, nreads=nreads.learn,
                                      errorEstimationFunction=dada2:::PacBioErrfun,
                                      multithread=multithread, BAND_SIZE=BAND_SIZE))
  err_plot <- plotErrors(err)
  pdf("plot_error_model.pdf", width=12, height=8, useDingbats=FALSE)
  print(err_plot)
  dev.off()

}else if(inp.dirR!='NULL'){#for paired read analysis

  err <- suppressWarnings(learnErrors(filts, nreads=nreads.learn, multithread=multithread))
  errR <- suppressWarnings(learnErrors(filtsR, nreads=nreads.learn, multithread=multithread))

}else{#for sinlge/pyro read analysis
  err <- suppressWarnings(learnErrors(filts, nreads=nreads.learn, multithread=multithread,
                                      HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY, BAND_SIZE=BAND_SIZE))
}

### PROCESS ALL SAMPLES ###
# Loop over rest in streaming fashion with learned error rates

if(inp.dirR =='NULL'){#for CCS/sinlge/pyro read analysis
  dds <- vector("list", length(filts))
  cat("4) Denoise samples ")
  cat("\n")
  for(j in seq(length(filts))) {
    drp <- derepFastq(filts[[j]])
    dds[[j]] <- dada(drp, err=err, multithread=multithread,HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY,
                     BAND_SIZE=BAND_SIZE, verbose=FALSE, OMEGA_C=1e-40)
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
                       priors=pseudo_priors, HOMOPOLYMER_GAP_PENALTY=HOMOPOLYMER_GAP_PENALTY,
                       BAND_SIZE=BAND_SIZE, verbose=FALSE, OMEGA_C=1e-40)
      cat(".")
    }
    cat("\n")
    ### \code copied from previous loop through samples in this script
  }
  # Make sequence table
  seqtab <- makeSequenceTable(dds)
}else{#for paired read analysis
  denoisedF <- rep(0, length(filts))
  ddsF <- vector("list", length(filts))
  ddsR <- vector("list", length(filts))
  mergers <- vector("list", length(filts))
  cat("3) Denoise samples ")

  for(j in seq(length(filts))) {
    drpF <- derepFastq(filts[[j]])
    ddsF[[j]] <- dada(drpF, err=err, multithread=multithread, verbose=FALSE, OMEGA_C=1e-40)
    drpR <- derepFastq(filtsR[[j]])
    ddsR[[j]] <- dada(drpR, err=errR, multithread=multithread, verbose=FALSE, OMEGA_C=1e-40)
    cat(".")
  }
  cat("\n")
  if(poolMethod == "pseudo") {
    cat("  Pseudo-pool step ")
    ### TEMPORARY, to be removed once 1.12 makes its way to Q2
    ### Needed for now to manage pseudo-pooling memory, as 1.10 didn't do this appropriately.
    ### pseudo_priors code copied from dada2.R
    stF <- makeSequenceTable(ddsF)
    pseudo_priorsF <- colnames(stF)[colSums(stF>0) >= 2 | colSums(stF) >= Inf]
    rm(stF)
    stR <- makeSequenceTable(ddsR)
    pseudo_priorsR <- colnames(stR)[colSums(stR>0) >= 2 | colSums(stR) >= Inf]
    rm(stR)
    ### \pseudo_priors code copied from dada2.R
    ### code copied from previous loop through samples in this script
    for(j in seq(length(filts))) {
      drpF <- derepFastq(filts[[j]])
      ddsF[[j]] <- dada(drpF, err=err, priors=pseudo_priorsF,
                        multithread=multithread, verbose=FALSE,OMEGA_C=1e-40)
      drpR <- derepFastq(filtsR[[j]])
      ddsR[[j]] <- dada(drpR, err=errR, priors=pseudo_priorsR,
                        multithread=multithread, verbose=FALSE,OMEGA_C=1e-40)
      cat(".")
    }
    cat("\n")
    ### \code copied from previous loop through samples in this script
  }

  ### Now loop through and do merging
  for(j in seq(length(filts))) {
    drpF <- derepFastq(filts[[j]])
    drpR <- derepFastq(filtsR[[j]])
    mergers[[j]] <- mergePairs(ddsF[[j]], drpF, ddsR[[j]], drpR, minOverlap=minOverlap)
    denoisedF[[j]] <- getN(ddsF[[j]])
    cat(".")
  }
  cat("\n")
  # Make sequence table
  seqtab <- makeSequenceTable(mergers)

}


### Remove chimeras
cat("5) Remove chimeras (method = ", chimeraMethod, ")\n", sep="")
if(chimeraMethod %in% c("pooled", "consensus")) {
  seqtab.nochim <- removeBimeraDenovo(seqtab, method=chimeraMethod, minFoldParentOverAbundance=minParentFold, allowOneOff=allowOneOff, multithread=multithread)
} else { # No chimera removal, copy seqtab to seqtab.nochim
  seqtab.nochim <- seqtab
}

### REPORT READ FRACTIONS THROUGH PIPELINE ###
cat("6) Report read numbers through the pipeline\n")
if(inp.dirR =='NULL'){
  # Handle edge cases: Samples lost in filtering; One sample
  if(primer.removed.dir!='NULL'){ #for CCS read analysis
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
  }else{#for paired end reads
    # Handle edge cases: Samples lost in filtering; One sample
    track <- cbind(out, matrix(0, nrow=nrow(out), ncol=3))
    colnames(track) <- c("input", "filtered", "denoised", "merged", "non-chimeric")
    passed.filtering <- track[,"filtered"] > 0
    track[passed.filtering,"denoised"] <- denoisedF
    track[passed.filtering,"merged"] <- rowSums(seqtab)
    track[passed.filtering,"non-chimeric"] <- rowSums(seqtab.nochim)
    write.table(track, out.track, sep="\t", row.names=TRUE, col.names=NA,
                quote=FALSE)
  }
}

### WRITE OUTPUT AND QUIT ###
# Formatting as tsv plain-text sequence table table
cat("7) Write output\n")
seqtab.nochim <- t(seqtab.nochim) # QIIME has OTUs as rows
col.names <- basename(filts)
col.names[[1]] <- paste0("#OTU ID\t", col.names[[1]])
write.table(seqtab.nochim, out.path, sep="\t",
            row.names=TRUE, col.names=col.names, quote=FALSE)
saveRDS(seqtab.nochim, "seqtab_nochim.rds") ### TESTING

q(status=0)
