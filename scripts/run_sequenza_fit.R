args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_sequenza_fit.R <binned.seqz.gz> <outdir> <sample_id>")
}
seqz_file <- args[[1]]
outdir <- args[[2]]
sample_id <- args[[3]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(sequenza))

patch_sequenza_gc_sample_stats <- function() {
  ns <- asNamespace("sequenza")
  current <- get("gc.sample.stats", envir = ns)
  body_text <- paste(deparse(body(current)), collapse = "\n")
  if (!grepl("parallel = parallel", body_text, fixed = TRUE)) {
    return(invisible(FALSE))
  }

  fixed <- function(file, col_types = "c--dd----d----", buffer = 33554432,
                    parallel = 2L, verbose = TRUE) {
    con <- gzfile(file, "rb")
    on.exit(close(con), add = TRUE)
    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, col_types) {
      x <- readr::read_tsv(
        file = paste(iotools::mstrsplit(x), collapse = "\n"),
        col_types = col_types,
        col_names = FALSE,
        skip = 0,
        n_max = Inf,
        progress = FALSE
      )
      u_chr <- unique(x[, 1])
      n_chr <- table(x[, 1])
      gc1 <- lapply(split(x[, 2], x[, 4]), table)
      gc2 <- lapply(split(x[, 3], x[, 4]), table)
      if (verbose) {
        message(".", appendLF = FALSE)
      }
      list(unique = u_chr, lines = n_chr, gc_nor = gc1, gc_tum = gc2)
    }
    if (verbose) {
      message("Collecting GC information ", appendLF = FALSE)
    }
    res <- iotools::chunk.apply(
      input = con,
      FUN = parse_chunck,
      col_types = col_types,
      CH.MAX.SIZE = buffer,
      CH.PARALLEL = parallel
    )
    if (verbose) {
      message(" done\n")
    }
    get("unfold_gc", envir = asNamespace("sequenza"))(res, stats = TRUE)
  }
  environment(fixed) <- ns
  unlockBinding("gc.sample.stats", ns)
  assign("gc.sample.stats", fixed, envir = ns)
  lockBinding("gc.sample.stats", ns)
  invisible(TRUE)
}

patched <- patch_sequenza_gc_sample_stats()
cat(sprintf("[%s] sequenza gc.sample.stats compatibility patch applied: %s\n", Sys.time(), patched))

cat(sprintf("[%s] sequenza.extract %s\n", Sys.time(), seqz_file))
seqz_data <- sequenza.extract(seqz_file, verbose = FALSE)
cat(sprintf("[%s] sequenza.fit %s\n", Sys.time(), sample_id))
CP <- sequenza.fit(seqz_data)
saveRDS(seqz_data, file.path(outdir, paste0(sample_id, ".sequenza_extract.rds")))
saveRDS(CP, file.path(outdir, paste0(sample_id, ".sequenza_fit.rds")))

cat(sprintf("[%s] sequenza.results %s\n", Sys.time(), sample_id))
results_status <- tryCatch({
  sequenza.results(sequenza.extract = seqz_data, cp.table = CP, sample.id = sample_id, out.dir = outdir)
  "ok"
}, error = function(e) {
  msg <- conditionMessage(e)
  cat(sprintf("[%s] WARNING sequenza.results failed after partial output: %s\n", Sys.time(), msg))
  writeLines(msg, file.path(outdir, paste0(sample_id, ".sequenza_results.warning.txt")))
  "warning"
})

make_cp_table <- function(CP) {
  if (is.data.frame(CP)) {
    return(CP)
  }
  if (is.list(CP) && all(c("ploidy", "cellularity", "lpp") %in% names(CP))) {
    grid <- expand.grid(ploidy = CP$ploidy, cellularity = CP$cellularity)
    grid$lpp <- as.vector(CP$lpp)
    grid$posterior <- grid$lpp
    grid$sequenza_results_status <- results_status
    return(grid)
  }
  data.frame(raw_class = paste(class(CP), collapse = ","), sequenza_results_status = results_status)
}

cp_table <- make_cp_table(CP)
cp_out <- file.path(outdir, paste0(sample_id, ".cp.table.tsv"))
write.table(cp_table, file = cp_out, sep = "\t", quote = FALSE, row.names = FALSE)

if ("posterior" %in% names(cp_table)) {
  best <- cp_table[which.max(cp_table$posterior), , drop = FALSE]
} else if ("score" %in% names(cp_table)) {
  best <- cp_table[order(cp_table$score), , drop = FALSE][1, , drop = FALSE]
} else {
  best <- cp_table[1, , drop = FALSE]
}
summary_out <- file.path(outdir, paste0(sample_id, ".sequenza_summary.tsv"))
write.table(best, file = summary_out, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("[%s] done %s\n", Sys.time(), sample_id))
