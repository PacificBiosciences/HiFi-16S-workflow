library(dada2)

# Read from command line argument the input FASTQ
args <- commandArgs(trailingOnly = TRUE)
fnFs <- args[1]
cpu <- args[2]

# Cap CPU to 16. Higher doesn't really help
if (cpu > 16) {
  cpu <- 16
}
# Learn the error rates
errF <- learnErrors(fnFs, multithread=cpu,  errorEstimationFunction=dada2:::PacBioErrfun)
err_plot <- plotErrors(errF)
pdf("plot_error_model.pdf", width=12, height=8, useDingbats=FALSE)
print(err_plot)
dev.off()

# Save as RDS
saveRDS(errF, file="errorfun.rds")
