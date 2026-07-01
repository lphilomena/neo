#!/usr/bin/env Rscript
# Fit FACETS from a snp-pileup CSV (normal then tumor BAM order).
# Usage: Rscript facets_fit_from_pileup.R <pileup.csv[.gz]> <output.rds>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: facets_fit_from_pileup.R <pileup.csv[.gz]> <output.rds>")
}

pileup_path <- args[1]
rds_path <- args[2]
ndepth <- as.integer(Sys.getenv("FACETS_NDEPTH", "5"))
cval_pre <- as.integer(Sys.getenv("FACETS_CVAL_PRE", "25"))
cval_proc <- as.integer(Sys.getenv("FACETS_CVAL_PROC", "25"))
min_nhet <- as.integer(Sys.getenv("FACETS_MIN_NHET", "5"))

suppressPackageStartupMessages(library(facets))

set.seed(1234)
rcmat <- readSnpMatrix(pileup_path)
cat("pileup SNPs:", nrow(rcmat), "\n")
if (nrow(rcmat) < 10) {
  stop("Too few SNPs in pileup (", nrow(rcmat), "); need deeper BAMs or lower -r threshold")
}

xx <- preProcSample(rcmat, gbuild = "hg38", ndepth = ndepth, ndepthmax = 2000, cval = cval_pre)
cat("preProcSample segments:", length(xx$chromlevels), "chromosomes\n")

oo <- procSample(xx, cval = cval_proc, min.nhet = min_nhet)
fit <- emcncf(oo, min.nhet = min_nhet)

saveRDS(fit, rds_path)
cat("purity:", fit$purity, "ploidy:", fit$ploidy, "\n")
if (!is.null(fit$emflags) && nzchar(fit$emflags)) {
  cat("emflags:", fit$emflags, "\n")
}
