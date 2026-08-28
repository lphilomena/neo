#!/usr/bin/env Rscript
# Export purity/ploidy and FACETS cncf segments from an emcncf RDS.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: facets_export_from_rds.R <facets.rds> <purity.tsv> <cncf.tsv>")
}

fit <- readRDS(args[1])
purity <- suppressWarnings(as.numeric(fit$purity))
ploidy <- suppressWarnings(as.numeric(fit$ploidy))
if (length(purity) != 1 || !is.finite(purity)) {
  stop("FACETS RDS does not contain a finite scalar purity")
}
if (length(ploidy) != 1 || !is.finite(ploidy)) {
  stop("FACETS RDS does not contain a finite scalar ploidy")
}

write.table(
  data.frame(purity = purity, ploidy = ploidy),
  file = args[2], sep = "\t", quote = FALSE, row.names = FALSE
)

cncf <- fit$cncf
if (is.null(cncf) || !is.data.frame(cncf) || nrow(cncf) == 0) {
  stop("FACETS RDS does not contain non-empty cncf segments")
}
write.table(cncf, file = args[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat("purity:", purity, "ploidy:", ploidy, "segments:", nrow(cncf), "\n")
