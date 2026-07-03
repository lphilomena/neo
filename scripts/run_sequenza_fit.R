args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_sequenza_fit.R <binned.seqz.gz> <outdir> <sample_id>")
}
seqz_file <- args[[1]]
outdir <- args[[2]]
sample_id <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(sequenza))

cat(sprintf("[%s] sequenza.extract %s\n", Sys.time(), seqz_file))
seqz_data <- sequenza.extract(seqz_file, verbose = FALSE)
cat(sprintf("[%s] sequenza.fit %s\n", Sys.time(), sample_id))
CP <- sequenza.fit(seqz_data)
saveRDS(seqz_data, file.path(outdir, paste0(sample_id, ".sequenza_extract.rds")))
saveRDS(CP, file.path(outdir, paste0(sample_id, ".sequenza_fit.rds")))

cat(sprintf("[%s] sequenza.results %s\n", Sys.time(), sample_id))
sequenza.results(sequenza.extract = seqz_data, cp.table = CP, sample.id = sample_id, out.dir = outdir)

cp_out <- file.path(outdir, paste0(sample_id, ".cp.table.tsv"))
write.table(CP, file = cp_out, sep = "\t", quote = FALSE, row.names = FALSE)

best <- CP[order(CP$score), , drop = FALSE][1, , drop = FALSE]
summary_out <- file.path(outdir, paste0(sample_id, ".sequenza_summary.tsv"))
write.table(best, file = summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[%s] done %s\n", Sys.time(), sample_id))
