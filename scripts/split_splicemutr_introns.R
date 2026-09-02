#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option("--input", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--chunks", type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = options))
if (is.null(opt$input) || is.null(opt$outdir) || opt$chunks < 1) {
  stop("--input, --outdir, and --chunks >= 1 are required")
}

introns <- readRDS(opt$input)
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
if (nrow(introns) == 0) {
  saveRDS(introns, file.path(opt$outdir, "introns.chunk_001.rds"))
  quit(status = 0)
}
chunk_count <- min(opt$chunks, nrow(introns))
groups <- split(seq_len(nrow(introns)), (seq_len(nrow(introns)) - 1) %% chunk_count)
for (index in seq_along(groups)) {
  saveRDS(introns[groups[[index]], , drop = FALSE],
          file.path(opt$outdir, sprintf("introns.chunk_%03d.rds", index)))
}
cat(sprintf("split %d introns into %d chunks\n", nrow(introns), chunk_count))
