library(dada2)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

seqs <- args[1]
threads <- as.numeric(args[2])
silva_db <- args[3]
gtdb_db <- args[4]
refseq_db <- args[5]
minBoot_num <- as.numeric(args[6])

seqs <- getSequences(seqs)
otu_id <- names(seqs)

# Assign with all three db
silva_spec <- assignTaxonomy(seqs, refFasta = silva_db, minBoot = minBoot_num,
                             multithread = threads, outputBootstraps = TRUE)
# Convert to df to handle n=1
silva_spec <- as.data.frame(silva_spec)
colnames(silva_spec)[grepl("tax.", colnames(silva_spec))] <-
  gsub("tax\\.(.*)", "\\1", colnames(silva_spec)[grepl("tax.", colnames(silva_spec))])
buf <- data.frame("Assignment" = rep("Silva 138.1", nrow(silva_spec)))
silva_spec <- cbind(silva_spec, buf)
write.table(data.frame("Feature ID" = otu_id, silva_spec), "silva_nb.tsv", quote = FALSE,
            sep = "\t",row.names = FALSE)

gtdb_spec <- assignTaxonomy(seqs, refFasta = gtdb_db, minBoot = minBoot_num,
                            multithread = threads, outputBootstraps = TRUE)
gtdb_spec <- as.data.frame(gtdb_spec)
colnames(gtdb_spec)[grepl("tax.", colnames(gtdb_spec))] <-
  gsub("tax\\.(.*)", "\\1", colnames(gtdb_spec)[grepl("tax.", colnames(gtdb_spec))])
# Rename species
gtdb_spec[, 'Species'] <- gsub("(.*)\\(.*\\)", "\\1", gtdb_spec[, 'Species'])
# Replace GTDB underscore with space
gtdb_spec[, 'Species'] <- gsub("_", "\\ ", gtdb_spec[, 'Species'])
buf <- data.frame("Assignment" = rep("GTDB r207", nrow(gtdb_spec)))
gtdb_spec <- cbind(gtdb_spec, buf)
write.table(data.frame("Feature ID" = otu_id, gtdb_spec), "gtdb_nb.tsv", quote = FALSE,
            sep = "\t",row.names = FALSE)

refseq_spec <- assignTaxonomy(seqs, refFasta = refseq_db, minBoot = minBoot_num,
                              multithread = threads, outputBootstraps = TRUE)
refseq_spec <- as.data.frame(refseq_spec)
colnames(refseq_spec)[grepl("tax.", colnames(refseq_spec))] <-
  gsub("tax\\.(.*)", "\\1", colnames(refseq_spec)[grepl("tax.", colnames(refseq_spec))])
# Rename species
refseq_spec[, 'Species'] <- gsub("(.*)_(.*?)_.*\\(.*\\)", "\\2", refseq_spec[, 'Species'])
buf <- data.frame("Assignment" = rep("RefSeq + RDP", nrow(refseq_spec)))
refseq_spec <- cbind(refseq_spec, buf)
write.table(data.frame("Feature ID" = otu_id, refseq_spec), "refseq_rdp_nb.tsv", quote = FALSE,
            sep = "\t",row.names = FALSE)

# Iteratively join at species level
final_spec <- gtdb_spec
final_spec[is.na(final_spec[, 'Species']), ] <- silva_spec[is.na(final_spec[, 'Species']), ]
final_spec[is.na(final_spec[, 'Species']), ] <- refseq_spec[is.na(final_spec[, 'Species']), ]
# For those that's been replaced by either SILVA or refSeq, paste species name
class_nonGTDB <- final_spec[is.na(gtdb_spec[, 'Species']), ]
class_nonGTDB <- class_nonGTDB[!is.na(class_nonGTDB[, 'Species']),]
final_spec[rownames(class_nonGTDB), 'Species'] <-
  paste(final_spec[rownames(class_nonGTDB), 'Genus', drop=TRUE],
        final_spec[rownames(class_nonGTDB), 'Species', drop=TRUE], sep=" ")

# Iteratively join at genus level
final_spec[is.na(final_spec[, 'Genus']), ] <- silva_spec[is.na(final_spec[, 'Genus']), ]
final_spec[is.na(final_spec[, 'Genus']), ] <- refseq_spec[is.na(final_spec[, 'Genus']), ]

# uncultured or metagenome
final_spec[grepl("uncultured", final_spec[, 'Species'], ignore.case = TRUE)] <- 
  refseq_spec[grepl("uncultured", final_spec[, 'Species'], ignore.case = TRUE)]
final_spec[grepl("metagenome", final_spec[, 'Species'], ignore.case = TRUE)] <- 
  refseq_spec[grepl("metagenome", final_spec[, 'Species'], ignore.case = TRUE)]

# uncultured or metagenome at genus level
final_spec[grepl("uncultured", final_spec[, 'Genus'], ignore.case = TRUE)] <- 
  gtdb_spec[grepl("uncultured", final_spec[, 'Genus'], ignore.case = TRUE)]
final_spec[grepl("metagenome", final_spec[, 'Genus'], ignore.case = TRUE)] <- 
  gtdb_spec[grepl("metagenome", final_spec[, 'Genus'], ignore.case = TRUE)]

# If still NA, revert back to GTDB
final_spec[is.na(final_spec[, 'Species']), ] <- gtdb_spec[is.na(final_spec[, 'Species']), ]
final_spec[is.na(final_spec[, 'Genus']), ] <- gtdb_spec[is.na(final_spec[, 'Genus']), ]
rownames(final_spec) <- otu_id

# For each row, look for the deepest level assigned and assign confidence value
to_save <- data.frame()
to_save2 <- data.frame()
for (i in 1:nrow(final_spec)){
  row <- final_spec[i, ]
  all_na <- which(is.na(row))
  if(length(all_na) != 0){
    min_ind <- min(all_na)
  } else {
    min_ind <- 14
  }
  # Max level 7. If no NA, set to species confidence (col 14)
  if (min_ind <= 7){
    conf <- row[[min_ind + 7]]
  } else {
    conf <- row[[min_ind]]
  }
  # If NA, set to unclassified
  for (j in all_na){
    if (j <= 7){
      row[j] <- "Unclassified"
    }
  }
  tosave_row <- data.frame(
    "Feature ID" = rownames(row),
    "Taxon" = paste(
      paste(c("d__", "p__", "c__", "o__", "f__", "g__", "s__"), c(row[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")]), sep=""), 
      collapse = "; "),
    "Confidence" = conf,
    "Assignment Database" = row['Assignment']
  )
  to_save <- rbind(to_save, tosave_row[, 1:3])
  to_save2 <- rbind(to_save2, tosave_row)
}

write.table(to_save, "best_taxonomy.tsv", quote = FALSE,
            sep = "\t",row.names = FALSE, col.names = c("Feature ID", "Taxon", "Confidence"))

write.table(to_save2, "best_taxonomy_withDB.tsv", quote = FALSE,
           sep = "\t",row.names = FALSE, 
           col.names = c("Feature ID", "Taxon", "Confidence", "Assignment Database"))
