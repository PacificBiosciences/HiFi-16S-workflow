library(dada2)

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) != 7) {
  stop("Wrong number of arguments", call. = FALSE)
}

seqs <- args[1]
threads <- as.numeric(args[2])
silva_db <- args[3]
gtdb_db <- args[4]
gg2_db <- args[5]
minBoot_num <- as.numeric(args[6])

# Allow specifying which database to prioritize
db_to_prioritize <- args[7]
# If not specified, prioritize GreenGenes2
if (is.na(db_to_prioritize)) {
  db_to_prioritize <- "GG2"
}
# Make sure db_to_prioritize is one of the three options
if (!(db_to_prioritize %in% c("GG2", "GTDB", "Silva"))) {
  stop("Invalid (Typo?) database to prioritize (Must be either GG2, GTDB or Silva)", call. = FALSE)
}

seqs <- getSequences(seqs)
otu_id <- names(seqs)

# Assign with all three db
silva_spec <- assignTaxonomy(seqs,
  refFasta = silva_db, minBoot = minBoot_num,
  multithread = threads, outputBootstraps = TRUE
)
# Convert to df to handle n=1
silva_spec <- as.data.frame(silva_spec)
colnames(silva_spec)[grepl("tax.", colnames(silva_spec))] <-
  gsub("tax\\.(.*)", "\\1", colnames(silva_spec)[grepl("tax.", colnames(silva_spec))])
buf <- data.frame("Assignment" = rep("Silva 138.1", nrow(silva_spec)))
silva_spec <- cbind(silva_spec, buf)
# Paste genus and species name to be species (Only if species is not NA)
silva_spec[!is.na(silva_spec[, "Species"]), "Species"] <-
  paste(silva_spec[!is.na(silva_spec[, "Species"]), "Genus"],
    silva_spec[!is.na(silva_spec[, "Species"]), "Species"],
    sep = " "
  )
write.table(data.frame("Feature ID" = otu_id, silva_spec), "silva_nb.tsv",
  quote = FALSE,
  sep = "\t", row.names = FALSE
)

gtdb_spec <- assignTaxonomy(seqs,
  refFasta = gtdb_db, minBoot = minBoot_num,
  multithread = threads, outputBootstraps = TRUE
)
gtdb_spec <- as.data.frame(gtdb_spec)
colnames(gtdb_spec)[grepl("tax.", colnames(gtdb_spec))] <-
  gsub("tax\\.(.*)", "\\1", colnames(gtdb_spec)[grepl("tax.", colnames(gtdb_spec))])
# Rename species (remove the part in parentheses, GTDB has the exact strain of the species)
gtdb_spec[, "Species"] <- gsub("(.*)\\(.*", "\\1", gtdb_spec[, "Species"])
# Replace GTDB underscore with space
gtdb_spec[, "Species"] <- gsub("_", "\\ ", gtdb_spec[, "Species"])
buf <- data.frame("Assignment" = rep("GTDB r207", nrow(gtdb_spec)))
gtdb_spec <- cbind(gtdb_spec, buf)
write.table(data.frame("Feature ID" = otu_id, gtdb_spec), "gtdb_nb.tsv",
  quote = FALSE,
  sep = "\t", row.names = FALSE
)

gg2_spec <- assignTaxonomy(seqs,
  refFasta = gg2_db, minBoot = minBoot_num,
  multithread = threads, outputBootstraps = TRUE
)
gg2_spec <- as.data.frame(gg2_spec)
tax_cols_gg2 <- colnames(gg2_spec)[grepl("tax.", colnames(gg2_spec))]
# Remove everything until first underscore in tax column names
for (col in tax_cols_gg2) {
  gg2_spec[[col]] <- gsub("^.*__(.*)", "\\1", gg2_spec[[col]])
}

colnames(gg2_spec)[grepl("tax.", colnames(gg2_spec))] <-
  gsub("tax\\.(.*)", "\\1", tax_cols_gg2)
# Rename species
buf <- data.frame("Assignment" = rep("GreenGenes2", nrow(gg2_spec)))
gg2_spec <- cbind(gg2_spec, buf)
# Paste genus and species name to be species (Only if species is not NA)
gg2_spec[!is.na(gg2_spec[, "Species"]), "Species"] <-
  paste(gg2_spec[!is.na(gg2_spec[, "Species"]), "Genus"],
    gg2_spec[!is.na(gg2_spec[, "Species"]), "Species"],
    sep = " "
  )

write.table(data.frame("Feature ID" = otu_id, gg2_spec), "gg2_nb.tsv",
  quote = FALSE,
  sep = "\t", row.names = FALSE
)

# Iteratively join at species level and Genus level
if (db_to_prioritize == "GTDB") {
  final_spec <- gtdb_spec
  final_spec[is.na(final_spec[, "Species"]), ] <-
    gg2_spec[is.na(final_spec[, "Species"]), ]
  final_spec[is.na(final_spec[, "Species"]), ] <-
    silva_spec[is.na(final_spec[, "Species"]), ]

  final_spec[is.na(final_spec[, "Genus"]), ] <-
    gtdb_spec[is.na(final_spec[, "Genus"]), ]
  final_spec[is.na(final_spec[, "Genus"]), ] <-
    gg2_spec[is.na(final_spec[, "Genus"]), ]
  final_spec[is.na(final_spec[, "Genus"]), ] <-
    silva_spec[is.na(final_spec[, "Genus"]), ]
} else if (db_to_prioritize == "GG2") {
  final_spec <- gg2_spec
  final_spec[is.na(final_spec[, "Species"]), ] <-
    gtdb_spec[is.na(final_spec[, "Species"]), ]
  final_spec[is.na(final_spec[, "Species"]), ] <-
    silva_spec[is.na(final_spec[, "Species"]), ]

  final_spec[is.na(final_spec[, "Genus"]), ] <-
    gtdb_spec[is.na(final_spec[, "Genus"]), ]
  final_spec[is.na(final_spec[, "Genus"]), ] <-
    silva_spec[is.na(final_spec[, "Genus"]), ]
} else if (db_to_prioritize == "Silva") {
  final_spec <- silva_spec
  final_spec[is.na(final_spec[, "Species"]), ] <-
    gg2_spec[is.na(final_spec[, "Species"]), ]
  final_spec[is.na(final_spec[, "Species"]), ] <-
    gtdb_spec[is.na(final_spec[, "Species"]), ]

  final_spec[is.na(final_spec[, "Genus"]), ] <-
    gg2_spec[is.na(final_spec[, "Genus"]), ]
  final_spec[is.na(final_spec[, "Genus"]), ] <-
    gtdb_spec[is.na(final_spec[, "Genus"]), ]
}

rownames(final_spec) <- otu_id

# For each row, look for the deepest level assigned and assign confidence value
to_save <- data.frame()
to_save2 <- data.frame()
for (i in 1:nrow(final_spec)) {
  row <- final_spec[i, ]
  all_na <- which(is.na(row))
  if (length(all_na) != 0) {
    min_ind <- min(all_na)
  } else {
    min_ind <- 14
  }
  # Max level 7. If no NA, set to species confidence (col 14)
  if (min_ind <= 7) {
    conf <- row[[min_ind + 7]]
  } else {
    conf <- row[[min_ind]]
  }
  # If NA, set to unclassified
  for (j in all_na) {
    if (j <= 7) {
      row[j] <- "Unclassified"
    }
  }
  tosave_row <- data.frame(
    "Feature ID" = rownames(row),
    "Taxon" = paste(
      paste(c("d__", "p__", "c__", "o__", "f__", "g__", "s__"),
        c(row[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]),
        sep = ""
      ),
      collapse = "; "
    ),
    "Confidence" = conf,
    "Assignment Database" = row["Assignment"]
  )
  to_save <- rbind(to_save, tosave_row[, 1:3])
  to_save2 <- rbind(to_save2, tosave_row)
}

write.table(to_save, "best_taxonomy.tsv",
  quote = FALSE,
  sep = "\t", row.names = FALSE, col.names = c("Feature ID", "Taxon", "Confidence")
)

write.table(to_save2, "best_taxonomy_withDB.tsv",
  quote = FALSE,
  sep = "\t", row.names = FALSE,
  col.names = c("Feature ID", "Taxon", "Confidence", "Assignment Database")
)
