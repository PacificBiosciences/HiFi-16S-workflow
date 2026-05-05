#!/usr/bin/env Rscript

# ---- helper: get script directory ----
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)

  if (length(file_arg) == 0) {
    return(getwd())
  }

  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}

# ---- load functions ----
script_dir <- get_script_dir()
source(file.path(script_dir, "lib", "taxonomy_harmonisation.R"))

# ---- args ----
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop(
    "Usage: harmonise_taxonomy.R <input.tsv> <db_name> <output.tsv>",
    call. = FALSE
  )
}

input_tsv  <- args[1]
db_name    <- tolower(args[2])
output_tsv <- args[3]

# ---- read ----
dt <- read.delim(
  input_tsv,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (!"Taxon" %in% names(dt)) {
  stop("Input file must contain a 'Taxon' column", call. = FALSE)
}

# ---- harmonise ----
if (db_name == "gg2") {
  dt$Taxon <- harmonise_gg2_taxonomy(dt$Taxon)
} else if (db_name == "gtdb") {
  dt$Taxon <- harmonise_gtdb_taxonomy(dt$Taxon)
} else if (db_name == "silva") {
  dt$Taxon <- harmonise_silva_taxonomy(dt$Taxon)
} else {
  message("Skipping taxonomy harmonisation for unsupported database: ", db_name)
}
# ---- write ----
write.table(
  dt,
  file = output_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)
