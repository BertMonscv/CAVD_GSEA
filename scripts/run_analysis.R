#!/usr/bin/env Rscript
# Execute the preserved analysis in a separate directory. No frozen file is edited.

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  usage <- paste(
    "Usage: Rscript --vanilla scripts/run_analysis.R --output-dir DIR",
    "       [--data-dir DIR] [--mode gsea|baseline|all] [--library DIR]",
    "The selected output subdirectory must be empty. No historical result is overwritten.",
    sep = "\n"
  )
  if (any(args %in% c("--help", "-h"))) { cat(usage, "\n"); return(invisible(NULL)) }
  cli <- commandArgs(trailingOnly = FALSE)
  script_arg <- cli[grepl("^--file=", cli)]
  if (length(script_arg) != 1L) stop("Run this wrapper with Rscript.")
  root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."), mustWork = TRUE)
  opts <- list(data = file.path(root, "csv"), output = NULL, mode = "gsea",
               library = file.path(root, "environment/library"))
  while (length(args)) {
    if (length(args) < 2L) stop(usage)
    name <- args[1]; value <- args[2]; args <- args[-c(1, 2)]
    if (name == "--data-dir") opts$data <- value
    else if (name == "--output-dir") opts$output <- value
    else if (name == "--mode") opts$mode <- value
    else if (name == "--library") opts$library <- value
    else stop("Unknown argument: ", name)
  }
  if (is.null(opts$output) || !opts$mode %in% c("gsea", "baseline", "all")) stop(usage)
  if (dir.exists(opts$library)) .libPaths(c(normalizePath(opts$library), .libPaths()))
  if (!requireNamespace("digest", quietly = TRUE)) stop("Install the recorded R environment first; digest is missing.")
  opts$data <- normalizePath(opts$data, mustWork = TRUE)
  dir.create(opts$output, recursive = TRUE, showWarnings = FALSE)
  opts$output <- normalizePath(opts$output, mustWork = TRUE)
  protected <- normalizePath(file.path(root, c("pre", "results")), mustWork = TRUE)
  if (any(vapply(protected, function(p) opts$output == p || startsWith(opts$output, paste0(p, "/")), logical(1)))) {
    stop("Choose an output directory outside pre/ and results/.")
  }
  manifest <- read.delim(file.path(root, "environment/geo_inputs.tsv"), stringsAsFactors = FALSE)
  hash <- function(path) digest::digest(file = path, algo = "sha256")
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(opts$data, manifest$filename[i])
    if (!file.exists(path) || hash(path) != manifest$sha256[i]) stop("Input checksum mismatch: ", path)
  }
  frozen <- read.delim(file.path(root, "environment/frozen_sources.tsv"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(frozen))) {
    if (hash(file.path(root, frozen$path[i])) != frozen$sha256[i]) stop("Frozen source checksum mismatch: ", frozen$path[i])
  }
  patch_assignments <- function(expressions, expected, replacement) {
    seen <- setNames(integer(length(expected)), names(expected))
    for (i in seq_along(expressions)) {
      x <- expressions[[i]]
      if (!is.call(x) || !identical(x[[1]], as.name("<-")) || !is.symbol(x[[2]])) next
      key <- as.character(x[[2]])
      if (!key %in% names(expected)) next
      if (!identical(x[[3]], expected[[key]])) stop("Unexpected source assignment for ", key)
      x[[3]] <- replacement[[key]]
      expressions[[i]] <- x
      seen[key] <- seen[key] + 1L
    }
    if (any(seen != 1L)) stop("Expected exactly one assignment for every wrapper override")
    expressions
  }
  run_one <- function(mode) {
    destination <- file.path(opts$output, mode)
    if (dir.exists(destination) && length(list.files(destination, all.files = TRUE, no.. = TRUE))) {
      stop("Output directory is not empty: ", destination)
    }
    dir.create(destination, recursive = TRUE, showWarnings = FALSE)
    path_keys <- c(".PATH_GSE51472", ".PATH_GSE12644", ".PATH_GSE83453")
    input_files <- paste0(c("GSE51472", "GSE12644", "GSE83453"), "_series_matrix.txt")
    expected <- as.list(setNames(file.path("csv", input_files), path_keys))
    replacement <- as.list(setNames(file.path(opts$data, input_files), path_keys))
    if (mode == "gsea") {
      source <- file.path(root, "pre/run_gsea_v2.R")
      expected <- c(expected, list(RUN_BOOTSTRAP = TRUE, .PATH_GENESETS = "pre/Supp_Table_S1_v2.tsv",
                                   .PATH_OUT_DIR = "results", .PATH_SESSIONINFO = "pre/sessionInfo_v2.txt"))
      replacement <- c(replacement, list(RUN_BOOTSTRAP = FALSE,
                                         .PATH_GENESETS = file.path(root, "pre/Supp_Table_S1_v2.tsv"),
                                         .PATH_OUT_DIR = destination,
                                         .PATH_SESSIONINFO = file.path(destination, "sessionInfo.txt")))
    } else {
      source <- file.path(root, "pre/draw_figure1.R")
      packages <- c("limma", "sva", "ggplot2", "AnnotationDbi", "hgu133plus2.db", "illuminaHumanv4.db",
                    "clusterProfiler", "org.Hs.eg.db", "patchwork", "ggrepel")
      absent <- packages[!vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
      if (length(absent)) stop("Install baseline packages first: ", paste(absent, collapse = ", "))
      expected <- c(expected, list(.OUT = "results/Figure1_panels"))
      replacement <- c(replacement, list(.OUT = destination))
    }
    expressions <- patch_assignments(parse(file = source, keep.source = FALSE), expected, replacement)
    # Record only paths and bootstrap changes. All analysis expressions remain unchanged.
    write.table(data.frame(assignment = names(expected),
                           frozen = vapply(expected, function(x) paste(deparse(x), collapse = " "), character(1)),
                           wrapper = vapply(replacement, function(x) paste(deparse(x), collapse = " "), character(1))),
                file.path(destination, "wrapper_overrides.tsv"), sep = "\t", row.names = FALSE, quote = TRUE)
    writeLines(c(paste("started_utc:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
                 paste("source:", basename(source)), paste("source_sha256:", hash(source)),
                 paste("mode:", mode), "source files are unchanged"), file.path(destination, "run_provenance.txt"))
    execution <- new.env(parent = globalenv())
    eval(expressions, envir = execution)
    if (mode == "gsea") {
      # The archived schema includes CI columns; the submission table omits them.
      result <- read.delim(file.path(destination, "Supp_Table_S2_v2.tsv"), check.names = FALSE)
      result <- result[, !names(result) %in% c("NES_SD", "NES_lower_95CI", "NES_upper_95CI"), drop = FALSE]
      write.table(result, file.path(destination, "GSEA_results_no_CI.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    }
    writeLines(capture.output(sessionInfo()), file.path(destination, "sessionInfo.txt"))
    cat("\nSUCCESS: ", mode, " outputs: ", destination, "\n", sep = "")
  }
  for (mode in if (opts$mode == "all") c("gsea", "baseline") else opts$mode) run_one(mode)
}

main()
