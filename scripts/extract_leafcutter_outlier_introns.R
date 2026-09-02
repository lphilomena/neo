#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))

options <- list(
  make_option("--pvalues", type = "character"),
  make_option("--sample", type = "character"),
  make_option("--threshold", type = "double", default = 0.05),
  make_option("--adjust-method", dest = "adjust_method", type = "character", default = "BH"),
  make_option("--out", type = "character")
)
opt <- parse_args(OptionParser(option_list = options))

if (is.null(opt$pvalues) || is.null(opt$sample) || is.null(opt$out)) {
  stop("--pvalues, --sample, and --out are required")
}

pvalues <- read.table(opt$pvalues, header = TRUE, check.names = FALSE,
                      sep = "\t", quote = "", comment.char = "")
if (!(opt$sample %in% colnames(pvalues))) {
  stop(sprintf("target sample %s is absent from %s", opt$sample, opt$pvalues))
}

scores <- pvalues[[opt$sample]]
if (tolower(opt$adjust_method) != "none") {
  scores <- p.adjust(scores, method = opt$adjust_method)
}
selected <- rownames(pvalues)[!is.na(scores) & scores <= opt$threshold]
parts <- strsplit(selected, ":", fixed = TRUE)
valid <- vapply(parts, length, integer(1)) >= 4
parts <- parts[valid]

if (length(parts) == 0) {
  introns <- data.frame(chr = character(), start = integer(), end = integer(),
                        strand = character(), stringsAsFactors = FALSE)
} else {
  introns <- unique(data.frame(
    chr = vapply(parts, `[[`, character(1), 1),
    start = as.integer(vapply(parts, `[[`, character(1), 2)),
    end = as.integer(vapply(parts, `[[`, character(1), 3)),
    strand = vapply(parts, function(x) tail(strsplit(x[[4]], "_", fixed = TRUE)[[1]], 1),
                    character(1)),
    stringsAsFactors = FALSE
  ))
}

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
saveRDS(introns, opt$out)
cat(sprintf("selected %d outlier junctions for %s at %s <= %g\n",
            nrow(introns), opt$sample, opt$adjust_method, opt$threshold))
