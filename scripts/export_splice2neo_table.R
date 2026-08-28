#!/usr/bin/env Rscript
# Export a splice2neo data.frame/tibble from RDS/RData to a stable TSV.
args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[[i + 1]]
}
input <- value("--input")
output <- value("--output")
object_name <- value("--object")
strict <- "--strict" %in% args
if (is.null(input) || is.null(output)) stop("Usage: export_splice2neo_table.R --input FILE --output TSV [--object NAME] [--strict]")
if (!file.exists(input)) stop(paste("Input not found:", input))
obj <- NULL
if (grepl("\\.[Rr][Dd][Ss]$", input)) {
  obj <- readRDS(input)
} else {
  env <- new.env(parent = emptyenv())
  loaded <- load(input, envir = env)
  if (!is.null(object_name)) {
    if (!exists(object_name, envir = env, inherits = FALSE)) stop(paste("Object not found:", object_name))
    obj <- get(object_name, envir = env)
  } else {
    data_frames <- loaded[vapply(loaded, function(x) is.data.frame(get(x, envir = env)), logical(1))]
    if (length(data_frames) != 1) stop("RData must contain exactly one data.frame/tibble unless --object is supplied")
    obj <- get(data_frames[[1]], envir = env)
  }
}
if (is.list(obj) && !is.data.frame(obj)) {
  candidates <- obj[vapply(obj, is.data.frame, logical(1))]
  if (length(candidates) != 1) stop("RDS list must contain exactly one data.frame/tibble")
  obj <- candidates[[1]]
}
if (!is.data.frame(obj)) stop("Selected object is not a data.frame/tibble")
required <- c("junc_id")
missing <- setdiff(required, colnames(obj))
if (length(missing) > 0 && strict) stop(paste("Missing required columns:", paste(missing, collapse = ",")))
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(obj, file = output, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat(output, "\n")
