#!/usr/bin/env Rscript
# Per-allele FilterSamReads step from LOHHLA mapping (count.events + Picard).
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("usage: lohhla_filter_alleles.R <regionDir> <bamId> <hlaBam> <gatkDir> <numMisMatch> <allele>...")
}
regionDir <- args[[1]]
bamId <- args[[2]]
hlaBam <- args[[3]]
gatkDir <- args[[4]]
numMisMatch <- as.numeric(args[[5]])
alleles <- args[6:length(args)]

suppressPackageStartupMessages(require(Rsamtools, quietly = TRUE))

samtools_bin <- Sys.getenv("SAMTOOLS_BIN", unset = "samtools")
java_bin <- Sys.getenv("JAVA_BIN", unset = "java")

samtools_modern <- function(bin) {
  lines <- suppressWarnings(
    system(paste(shQuote(bin), "--version 2>/dev/null"), intern = TRUE)
  )
  length(lines) > 0 && grepl("^samtools [0-9]", lines[[1]])
}

run_cmd <- function(cmd, label) {
  status <- system(cmd)
  if (!is.na(status) && status != 0) {
    stop(label, " failed (exit ", status, "): ", cmd)
  }
  invisible(status)
}

picard.run <- function(tool, args) {
  jar <- file.path(gatkDir, "picard.jar")
  jar_arg <- if (file.exists(jar)) {
    paste0(shQuote(jar))
  } else {
    paste0(shQuote(file.path(gatkDir, paste0(tool, ".jar"))))
  }
  paste(java_bin, "-jar", jar_arg, tool, args)
}

count.events <- function(BAMfile, n) {
  x <- scanBam(BAMfile, index = BAMfile, param = ScanBamParam(what = scanBamWhat(), tag = "NM"))
  readIDs <- x[[1]][["qname"]]
  if (length(readIDs) == 0) {
    return(character())
  }
  cigar <- as.character(x[[1]][["cigar"]])
  editDistance <- x[[1]][["tag"]][["NM"]]
  if (is.null(editDistance)) {
    editDistance <- rep(0L, length(readIDs))
  } else {
    editDistance <- unlist(editDistance)
    if (length(editDistance) == 0) {
      editDistance <- rep(0L, length(readIDs))
    }
  }
  insertionCount <- sapply(cigar, function(boop) {
    length(grep(pattern = "I", x = unlist(strsplit(boop, split = ""))))
  })
  deletionCount <- sapply(cigar, function(boop) {
    length(grep(pattern = "D", x = unlist(strsplit(boop, split = ""))))
  })
  indelTotals <- sapply(cigar, function(boop) {
    tmp <- unlist(strsplit(gsub("([0-9]+)", "~\\1~", boop), "~"))
    Is <- grep(pattern = "I", x = tmp)
    Ds <- grep(pattern = "D", x = tmp)
    sum(as.numeric(tmp[(Is - 1)])) + sum(as.numeric(tmp[Ds - 1]))
  })
  misMatchCount <- editDistance - indelTotals
  eventCount <- misMatchCount + insertionCount + deletionCount
  names(eventCount) <- seq_along(eventCount)
  passed <- eventCount[which(eventCount <= n)]
  y <- readIDs[as.numeric(names(passed))]
  y <- names(table(y)[which(table(y) == 2)])
  y
}

for (allele in alleles) {
  outBam <- file.path(regionDir, paste0(bamId, ".type.", allele, ".filtered.bam"))
  if (file.exists(outBam) && file.info(outBam)$size > 0) {
    message("skip existing ", outBam)
    next
  }
  message("allele filter: ", allele)
  tempBam <- file.path(regionDir, paste0(bamId, ".temp.", allele, ".bam"))
  typeBam <- file.path(regionDir, paste0(bamId, ".type.", allele, ".bam"))
  typePrefix <- sub("\\.bam$", "", typeBam)
  passed.reads.file <- file.path(regionDir, paste0(bamId, ".", allele, ".passed.reads.txt"))

  # samtools 0.1.x (polysolver) lacks "sort -o"; redirect view output for portability.
  cmd <- paste(
    shQuote(samtools_bin), "view -b", shQuote(hlaBam), allele,
    ">", shQuote(tempBam)
  )
  run_cmd(cmd, "samtools view")
  if (samtools_modern(samtools_bin)) {
    cmd <- paste(
      shQuote(samtools_bin), "sort -o", shQuote(typeBam), shQuote(tempBam)
    )
  } else {
    cmd <- paste(shQuote(samtools_bin), "sort", shQuote(tempBam), shQuote(typePrefix))
  }
  run_cmd(cmd, "samtools sort")
  if (!file.exists(typeBam) || file.info(typeBam)$size == 0) {
    stop("expected sorted BAM missing: ", typeBam)
  }
  cmd <- paste(shQuote(samtools_bin), "index", shQuote(typeBam))
  run_cmd(cmd, "samtools index")

  passed.reads <- count.events(typeBam, n = numMisMatch)
  if (length(passed.reads) == 0) {
    file.copy(typeBam, outBam, overwrite = TRUE)
  } else {
    write.table(
      passed.reads,
      file = passed.reads.file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = FALSE
    )
    extractCMD <- picard.run(
      "FilterSamReads",
      paste0(
        " I=", typeBam,
        " FILTER=includeReadList READ_LIST_FILE=", passed.reads.file,
        " OUTPUT=", outBam
      )
    )
    run_cmd(extractCMD, "FilterSamReads")
  }
  cmd <- paste(shQuote(samtools_bin), "index", shQuote(outBam))
  run_cmd(cmd, "samtools index filtered BAM")

  unlink(c(tempBam, paste0(tempBam, ".bai"), typeBam, paste0(typeBam, ".bai")))
}
